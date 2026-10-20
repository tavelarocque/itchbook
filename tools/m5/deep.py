"""DeepLOB (Zhang, Zohren & Roberts 2019) and an order-flow variant, in PyTorch.

DeepLOB reads a (T x 40) window of the book: three convolution blocks that
first merge each level's price and size, then the two sides, then all ten
levels; an inception module over time; an LSTM; a classifier. The order-flow
variant (in the spirit of Kolm, Turiel & Westray 2023) reads (T x 12) windows
of multi-level OFI, signed volume and event counts, with the level-merging
convolutions replaced by one convolution across all inputs.

Both have one 3-class head (down, flat, up) per forecast horizon, trained
jointly, so one network serves all horizons.
"""

from __future__ import annotations

import torch
from torch import nn


def _conv(cin: int, cout: int, k: tuple[int, int], s: tuple[int, int] = (1, 1), pad: str | int = 0) -> list[nn.Module]:
    return [nn.Conv2d(cin, cout, k, stride=s, padding=pad), nn.LeakyReLU(0.01), nn.BatchNorm2d(cout)]


def _temporal(c: int) -> list[nn.Module]:
    # Two (4 x 1) convolutions over time. Each is padded with 1 step before and 2 after,
    # so the time axis keeps its length (Keras "same" for an even kernel).
    pad = [nn.ZeroPad2d((0, 0, 1, 2))]
    return pad + _conv(c, c, (4, 1)) + pad + _conv(c, c, (4, 1))


class Inception(nn.Module):
    def __init__(self, cin: int, c: int = 64) -> None:
        super().__init__()
        self.b1 = nn.Sequential(*_conv(cin, c, (1, 1), pad="same"), *_conv(c, c, (3, 1), pad="same"))
        self.b2 = nn.Sequential(*_conv(cin, c, (1, 1), pad="same"), *_conv(c, c, (5, 1), pad="same"))
        self.b3 = nn.Sequential(nn.MaxPool2d((3, 1), stride=(1, 1), padding=(1, 0)), *_conv(cin, c, (1, 1), pad="same"))

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return torch.cat([self.b1(x), self.b2(x), self.b3(x)], dim=1)


class DeepLOB(nn.Module):
    def __init__(self, n_features: int, n_heads: int, lob_layout: bool = True, filters: int = 32) -> None:
        super().__init__()
        f = filters
        if lob_layout:
            if n_features != 40:
                raise ValueError("the LOB layout expects 40 features")
            self.convs = nn.Sequential(
                *_conv(1, f, (1, 2), (1, 2)),  # price with size: 40 -> 20
                *_temporal(f),
                *_conv(f, f, (1, 2), (1, 2)),  # ask with bid: 20 -> 10
                *_temporal(f),
                *_conv(f, f, (1, 10)),  # all levels: 10 -> 1
                *_temporal(f),
            )
        else:
            self.convs = nn.Sequential(*_conv(1, f, (1, n_features)), *_temporal(f), *_temporal(f))
        self.inception = Inception(f)
        self.lstm = nn.LSTM(input_size=3 * 64, hidden_size=64, batch_first=True)
        self.heads = nn.ModuleList([nn.Linear(64, 3) for _ in range(n_heads)])

    def forward(self, x: torch.Tensor) -> list[torch.Tensor]:
        # x: (batch, T, features)
        z = self.convs(x.unsqueeze(1))  # (batch, f, T, 1)
        z = self.inception(z)  # (batch, 192, T, 1)
        z = z.squeeze(3).transpose(1, 2)  # (batch, T, 192)
        out, _ = self.lstm(z)
        h = out[:, -1, :]
        return [head(h) for head in self.heads]


def device() -> torch.device:
    if torch.backends.mps.is_available():
        return torch.device("mps")
    if torch.cuda.is_available():
        return torch.device("cuda")
    return torch.device("cpu")
