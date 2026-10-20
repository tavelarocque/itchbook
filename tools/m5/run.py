#!/usr/bin/env python3
"""M5: do deep order-book models beat an order-flow-imbalance baseline after costs?

  run.py universe --rank RANK.tsv --out study.json      fix the symbols (from the first day only)
  run.py ridge    --study study.json --data DIR --preds DIR
  run.py deep     --study study.json --data DIR --preds DIR --input lob|ofi [--epochs N --steps N]
  run.py evaluate --study study.json --data DIR --preds DIR --repo REPO

Days are split by date: models are fitted on the training days, every choice
(ridge penalty, decision thresholds, early stopping) is made on the
validation day, and the test days are touched once, by `evaluate`.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lob import HORIZONS, LOOKBACK, Study, SymbolDay  # noqa: E402

H_NAMES = list(HORIZONS)
H_STEPS = [HORIZONS[h] for h in H_NAMES]
TAKER_FEE_USD = 0.0030  # per share per side, NASDAQ standard fee to remove liquidity (2019-20)


def log(*a: object) -> None:
    print(time.strftime("%H:%M:%S"), *a, flush=True)


# ---- universe -------------------------------------------------------------------


def cmd_universe(a: argparse.Namespace) -> int:
    rows = [line.rstrip("\n").split("\t") for line in Path(a.rank).read_text().splitlines()]
    head, body = rows[0], rows[1:]
    i_sym, i_vwap = head.index("symbol"), head.index("vwap")
    picked = [r[i_sym] for r in body if float(r[i_vwap]) >= a.min_price][: a.n]
    study = {
        "rule": f"top {a.n} by trade prints on the first day, VWAP >= ${a.min_price:g}",
        "universe": picked,
        "train": a.train,
        "val": a.val,
        "test": a.test,
    }
    Path(a.out).write_text(json.dumps(study, indent=2) + "\n")
    print(json.dumps(study, indent=2))
    return 0


# ---- metrics --------------------------------------------------------------------


def macro_f1(true: np.ndarray, pred: np.ndarray) -> float:
    f = []
    for c in (0, 1, 2):
        tp = np.sum((pred == c) & (true == c))
        fp = np.sum((pred == c) & (true != c))
        fn = np.sum((pred != c) & (true == c))
        f.append(0.0 if tp == 0 else 2 * tp / (2 * tp + fp + fn))
    return float(np.mean(f))


def to_class(y: np.ndarray) -> np.ndarray:
    return (np.sign(y) + 1).astype(np.int8)  # 0 down, 1 flat, 2 up


# ---- ridge ----------------------------------------------------------------------


def cmd_ridge(a: argparse.Namespace) -> int:
    st = Study.from_json(Path(a.study), Path(a.data))
    scales = st.scales()
    k = None
    A = b = n = None
    for day in st.train:
        for s in st.universe:
            sd = st.load(day, s)
            X = sd.linear(scales[s])
            if k is None:
                k = X.shape[1] + 1
                A = np.zeros((len(H_STEPS), k, k))
                b = np.zeros((len(H_STEPS), k))
                n = np.zeros(len(H_STEPS))
            for j, h in enumerate(H_STEPS):
                y, ok = sd.targets(h)
                Z = np.hstack([np.ones((ok.sum(), 1)), X[ok]])
                A[j] += Z.T @ Z
                b[j] += Z.T @ y[ok]
                n[j] += ok.sum()
        log("ridge: accumulated", day)
    alphas = [0.0, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0]
    fits = {}
    for j, hn in enumerate(H_NAMES):
        mean = A[j][0, 1:] / n[j]
        cov = A[j][1:, 1:] / n[j] - np.outer(mean, mean)
        sd = np.sqrt(np.maximum(np.diag(cov), 1e-12))
        corr = cov / np.outer(sd, sd)
        ymean = b[j][0] / n[j]
        xy = (b[j][1:] / n[j] - mean * ymean) / sd
        fits[hn] = []
        for al in alphas:
            beta_std = np.linalg.solve(corr + al * np.eye(len(sd)), xy)
            beta = beta_std / sd
            fits[hn].append((al, beta, ymean - mean @ beta))
    # Choose the penalty on the validation day by R^2.
    chosen = {}
    for j, hn in enumerate(H_NAMES):
        best = None
        for al, beta, c in fits[hn]:
            sse = sst = 0.0
            for day in st.val:
                for s in st.universe:
                    sd_ = st.load(day, s)
                    y, ok = sd_.targets(H_STEPS[j])
                    p = sd_.linear(scales[s])[ok] @ beta + c
                    sse += float(np.sum((y[ok] - p) ** 2))
                    sst += float(np.sum((y[ok] - y[ok].mean()) ** 2))
            r2 = 1 - sse / sst
            if best is None or r2 > best[0]:
                best = (r2, al, beta, c)
        chosen[hn] = best
        log(f"ridge {hn}: alpha {best[1]:g}, validation R^2 {best[0]:.4f}")
    preds = Path(a.preds)
    preds.mkdir(parents=True, exist_ok=True)
    for day in st.days():
        for s in st.universe:
            sd_ = st.load(day, s)
            X = sd_.linear(scales[s])
            out = np.column_stack([X @ chosen[hn][2] + chosen[hn][3] for hn in H_NAMES]).astype(np.float32)
            np.save(preds / f"ridge.{day}.{s}.npy", out)
    meta = {hn: {"alpha": chosen[hn][1], "val_r2": chosen[hn][0], "coef": chosen[hn][2].tolist()} for hn in H_NAMES}
    (preds / "ridge.meta.json").write_text(json.dumps(meta, indent=1))
    return 0


# ---- deep -----------------------------------------------------------------------


class Windows:
    """All symbol-days of a split, normalised, with valid window end points."""

    def __init__(self, st: Study, days: list[str], view: str, scales: dict[str, float], stats=None) -> None:
        self.X, self.Y, self.ends, self.keys = [], [], [], []
        for day in days:
            for s in st.universe:
                sd = st.load(day, s)
                x = sd.lob() if view == "lob" else sd.ofi(scales[s])
                ys = [sd.targets(h) for h in H_STEPS]
                ok = np.logical_and.reduce([m for _, m in ys])
                ok[: LOOKBACK - 1] = False
                y = np.stack([to_class(t) for t, _ in ys], axis=1)
                self.X.append(x)
                self.Y.append(y)
                self.ends.append(np.flatnonzero(ok))
                self.keys.append((day, s))
        if stats is None:
            allx = np.concatenate([x[e] for x, e in zip(self.X, self.ends, strict=True)])
            stats = (allx.mean(axis=0), allx.std(axis=0) + 1e-6)
        self.stats = stats
        self.X = [((x - stats[0]) / stats[1]).astype(np.float16) for x in self.X]
        self.sizes = np.array([len(e) for e in self.ends])

    def sample(self, rng: np.random.Generator, n: int) -> tuple[np.ndarray, np.ndarray]:
        which = rng.choice(len(self.X), size=n, p=self.sizes / self.sizes.sum())
        xb = np.empty((n, LOOKBACK, self.X[0].shape[1]), dtype=np.float32)
        yb = np.empty((n, len(H_STEPS)), dtype=np.int64)
        for i in np.unique(which):
            sel = np.flatnonzero(which == i)
            t = self.ends[i][rng.integers(0, len(self.ends[i]), size=len(sel))]
            idx = t[:, None] + np.arange(-LOOKBACK + 1, 1)[None, :]
            xb[sel] = self.X[i][idx]
            yb[sel] = self.Y[i][t]
        return xb, yb


def cmd_deep(a: argparse.Namespace) -> int:
    import torch

    from deep import DeepLOB, device

    st = Study.from_json(Path(a.study), Path(a.data))
    scales = st.scales()
    dev = device()
    torch.manual_seed(a.seed)
    rng = np.random.default_rng(a.seed)
    log("deep: loading training days on", dev)
    tr = Windows(st, st.train, a.input, scales)
    va = Windows(st, st.val, a.input, scales, tr.stats)
    log(f"deep: {tr.sizes.sum():,} training windows, {va.sizes.sum():,} validation windows")
    nf = tr.X[0].shape[1]
    model = DeepLOB(nf, len(H_STEPS), lob_layout=a.input == "lob").to(dev)
    opt = torch.optim.Adam(model.parameters(), lr=a.lr)
    loss_fn = torch.nn.CrossEntropyLoss()
    vx, vy = va.sample(np.random.default_rng(12345), a.val_samples)
    name = "deeplob" if a.input == "lob" else "deepofi"
    preds = Path(a.preds)
    preds.mkdir(parents=True, exist_ok=True)
    ckpt = preds / f"{name}.pt"

    def predict(x: np.ndarray) -> np.ndarray:
        model.eval()
        out = []
        with torch.no_grad():
            for i in range(0, len(x), 4096):
                logits = model(torch.from_numpy(x[i : i + 4096]).to(dev))
                out.append(np.stack([torch.softmax(z, 1).cpu().numpy() for z in logits], axis=1))
        return np.concatenate(out)

    best, history, bad = -1.0, [], 0
    for epoch in range(a.epochs):
        model.train()
        t0 = time.time()
        tot = 0.0
        for step in range(a.steps):
            xb, yb = tr.sample(rng, a.batch)
            xb_t = torch.from_numpy(xb).to(dev)
            yb_t = torch.from_numpy(yb).to(dev)
            out = model(xb_t)
            loss = sum(loss_fn(o, yb_t[:, j]) for j, o in enumerate(out))
            opt.zero_grad()
            loss.backward()
            opt.step()
            tot += float(loss.detach())
            if step % 500 == 0:
                rate = (step + 1) * a.batch / (time.time() - t0)
                log(f"  epoch {epoch} step {step} loss {tot / (step + 1):.4f} ({rate:.0f} windows/s)")
        p = predict(vx)
        f1 = [macro_f1(vy[:, j], p[:, j].argmax(1)) for j in range(len(H_STEPS))]
        score = float(np.mean(f1))
        history.append({"epoch": epoch, "train_loss": tot / a.steps, "val_macro_f1": f1, "seconds": time.time() - t0})
        log(f"deep {name} epoch {epoch}: loss {tot / a.steps:.4f}, validation macro-F1 {np.round(f1, 4)}")
        if score > best:
            best, bad = score, 0
            torch.save(model.state_dict(), ckpt)
        else:
            bad += 1
            if bad >= a.patience:
                log("deep: early stop")
                break
    model.load_state_dict(torch.load(ckpt, map_location=dev))
    log("deep: predicting every day with the best epoch")
    for day in st.days():
        ws = Windows(st, [day], a.input, scales, tr.stats)
        for (d, s), x in zip(ws.keys, ws.X, strict=True):
            n = len(x)
            ends = np.arange(LOOKBACK - 1, n)
            out = np.full((n, len(H_STEPS), 3), np.nan, dtype=np.float16)
            for i in range(0, len(ends), 16384):
                e = ends[i : i + 16384]
                xb = x[e[:, None] + np.arange(-LOOKBACK + 1, 1)[None, :]].astype(np.float32)
                out[e] = predict(xb).astype(np.float16)
            np.save(preds / f"{name}.{d}.{s}.npy", out)
        log("  predicted", day)
    meta = {"input": a.input, "history": history, "best_val_macro_f1": best, "args": vars(a)}
    (preds / f"{name}.meta.json").write_text(json.dumps(meta, indent=1, default=str))
    return 0


# ---- evaluation -----------------------------------------------------------------


def scores(preds: Path, model: str, day: str, sym: str) -> np.ndarray:
    """(rows, horizons) score: expected direction, in ticks for ridge, P(up) - P(down) for deep models."""
    p = np.load(preds / f"{model}.{day}.{sym}.npy").astype(np.float32)
    return p if p.ndim == 2 else p[:, :, 2] - p[:, :, 0]


def classes(preds: Path, model: str, day: str, sym: str, dead: list[float]) -> np.ndarray:
    p = np.load(preds / f"{model}.{day}.{sym}.npy").astype(np.float32)
    if p.ndim == 2:  # ridge: a dead zone around zero is the flat class
        return np.stack(
            [np.where(p[:, j] > dead[j], 2, np.where(p[:, j] < -dead[j], 0, 1)) for j in range(p.shape[1])], 1
        )
    return np.nan_to_num(p, nan=0).argmax(2)


def trades(sd: SymbolDay, s: np.ndarray, h: int, theta: float, fee_ticks: np.ndarray) -> dict[str, np.ndarray]:
    """Taker round trips at decision points spaced h apart: buy at the ask if s > theta, sell at the bid
    if s < -theta, unwind at the opposite touch h steps later. Returns per-trade components in ticks."""
    n = len(sd.raw)
    v = sd.valid()
    t = np.arange(LOOKBACK, n - h, h)
    t = t[v[t] & v[t + h] & ~np.isnan(s[t])]
    side = np.where(s[t] > theta, 1, np.where(s[t] < -theta, -1, 0))
    t, side = t[side != 0], side[side != 0]
    tick = sd.tick()[t]
    bid, ask, mid = sd.bid, sd.ask, sd.mid()
    gross = side * (mid[t + h] - mid[t]) / tick
    spread_cost = ((ask[t] - bid[t]) / 2 + (ask[t + h] - bid[t + h]) / 2) / tick
    fees = 2 * fee_ticks[t]
    return {"gross": gross, "spread": spread_cost, "fees": fees, "net": gross - spread_cost - fees}


def cmd_evaluate(a: argparse.Namespace) -> int:
    import pandas as pd

    st = Study.from_json(Path(a.study), Path(a.data))
    preds = Path(a.preds)
    res = Path(a.repo) / "results"
    res.mkdir(exist_ok=True)
    models = [m for m in ("ridge", "deepofi", "deeplob") if (preds / f"{m}.{st.test[0]}.{st.universe[0]}.npy").exists()]
    cache: dict[tuple[str, str], SymbolDay] = {}

    def sd_of(day: str, s: str) -> SymbolDay:
        if (day, s) not in cache:
            cache[(day, s)] = st.load(day, s)
        return cache[(day, s)]

    # Ridge dead zones for classification, chosen on validation by macro-F1.
    dead = [0.0] * len(H_STEPS)
    if "ridge" in models:
        for j, h in enumerate(H_STEPS):
            best = (-1, 0.0)
            for d in [0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5]:
                tr_, pr_ = [], []
                for day in st.val:
                    for s in st.universe:
                        y, ok = sd_of(day, s).targets(h)
                        p = scores(preds, "ridge", day, s)[:, j]
                        tr_.append(to_class(y[ok]))
                        pr_.append(np.where(p[ok] > d, 2, np.where(p[ok] < -d, 0, 1)))
                f1 = macro_f1(np.concatenate(tr_), np.concatenate(pr_))
                best = max(best, (f1, d))
            dead[j] = best[1]

    rows = []
    for m in models + ["majority"]:
        for j, (hn, h) in enumerate(zip(H_NAMES, H_STEPS, strict=True)):
            tr_, pr_ = [], []
            for day in st.test:
                for s in st.universe:
                    y, ok = sd_of(day, s).targets(h)
                    ok = ok.copy()
                    ok[:LOOKBACK] = False
                    tr_.append(to_class(y[ok]))
                    if m == "majority":
                        pr_.append(np.ones(ok.sum(), dtype=np.int8))
                    else:
                        pr_.append(classes(preds, m, day, s, dead)[ok, j])
            t_, p_ = np.concatenate(tr_), np.concatenate(pr_)
            moved = t_ != 1
            rows.append(
                {
                    "model": m,
                    "horizon": hn,
                    "accuracy": float(np.mean(t_ == p_)),
                    "macro_f1": macro_f1(t_, p_),
                    "direction_hit_when_moved_and_called": float(
                        np.mean(t_[moved & (p_ != 1)] == p_[moved & (p_ != 1)])
                    )
                    if np.any(moved & (p_ != 1))
                    else np.nan,
                    "share_flat_true": float(np.mean(t_ == 1)),
                }
            )
    cls = pd.DataFrame(rows)
    cls.to_csv(res / "m5_classification.tsv", sep="\t", index=False, float_format="%.4f")
    print(cls.to_string(index=False))

    # Trading: thresholds chosen on validation to maximise total net P&L, then applied once to test.
    def fee_ticks(sd: SymbolDay, fee: float) -> np.ndarray:
        return fee / (sd.tick() / 1e4)

    trows = []
    for fee_name, fee in [("taker fee", TAKER_FEE_USD), ("no fee", 0.0)]:
        for m in models + ["oracle"]:
            for j, (hn, h) in enumerate(zip(H_NAMES, H_STEPS, strict=True)):

                def score_of(day: str, s: str, m=m, j=j, h=h) -> np.ndarray:
                    if m == "oracle":
                        y, _ = sd_of(day, s).targets(h)
                        return y
                    return scores(preds, m, day, s)[:, j]

                grid = np.unique(
                    np.quantile(
                        np.abs(np.concatenate([score_of(d, s)[LOOKBACK:] for d in st.val for s in st.universe])),
                        [0.0, 0.5, 0.8, 0.9, 0.95, 0.98, 0.99, 0.995, 0.999],
                    )
                )
                best = (0.0, None)  # no trading is always available and earns zero
                for th in grid:
                    tot = 0.0
                    for day in st.val:
                        for s in st.universe:
                            sd = sd_of(day, s)
                            tot += float(trades(sd, score_of(day, s), h, th, fee_ticks(sd, fee))["net"].sum())
                    if tot > best[0]:
                        best = (tot, th)
                parts = {k: [] for k in ("gross", "spread", "fees", "net")}
                if best[1] is not None:
                    for day in st.test:
                        for s in st.universe:
                            sd = sd_of(day, s)
                            tr = trades(sd, score_of(day, s), h, best[1], fee_ticks(sd, fee))
                            for k in parts:
                                parts[k].append(tr[k])
                cat = {k: np.concatenate(v) if v else np.zeros(0) for k, v in parts.items()}
                trows.append(
                    {
                        "costs": fee_name,
                        "model": m,
                        "horizon": hn,
                        "val_net_ticks": best[0],
                        "threshold": best[1] if best[1] is not None else np.nan,
                        "test_trades": len(cat["net"]),
                        "gross_per_trade": cat["gross"].mean() if len(cat["net"]) else np.nan,
                        "spread_per_trade": cat["spread"].mean() if len(cat["net"]) else np.nan,
                        "fees_per_trade": cat["fees"].mean() if len(cat["net"]) else np.nan,
                        "net_per_trade": cat["net"].mean() if len(cat["net"]) else np.nan,
                        "test_net_ticks": cat["net"].sum(),
                        "hit_rate": float(np.mean(cat["gross"] > 0)) if len(cat["net"]) else np.nan,
                    }
                )
    tr = pd.DataFrame(trows)
    tr.to_csv(res / "m5_trading.tsv", sep="\t", index=False, float_format="%.4f")
    print(tr.to_string(index=False))
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    u = sub.add_parser("universe")
    u.add_argument("--rank", required=True)
    u.add_argument("--out", required=True)
    u.add_argument("--n", type=int, default=20)
    u.add_argument("--min-price", type=float, default=5.0)
    u.add_argument("--train", nargs="+", default=["01302019", "03272019", "07302019", "08302019"])
    u.add_argument("--val", nargs="+", default=["10302019"])
    u.add_argument("--test", nargs="+", default=["12302019", "01302020"])
    for name in ("ridge", "deep", "evaluate"):
        p = sub.add_parser(name)
        p.add_argument("--study", required=True)
        p.add_argument("--data", required=True)
        p.add_argument("--preds", required=True)
        if name == "deep":
            p.add_argument("--input", choices=["lob", "ofi"], required=True)
            p.add_argument("--epochs", type=int, default=8)
            p.add_argument("--steps", type=int, default=2000)
            p.add_argument("--batch", type=int, default=256)
            p.add_argument("--lr", type=float, default=1e-3)
            p.add_argument("--patience", type=int, default=2)
            p.add_argument("--val-samples", type=int, default=100_000)
            p.add_argument("--seed", type=int, default=1)
        if name == "evaluate":
            p.add_argument("--repo", default=".")
    a = ap.parse_args()
    return {"universe": cmd_universe, "ridge": cmd_ridge, "deep": cmd_deep, "evaluate": cmd_evaluate}[a.cmd](a)


if __name__ == "__main__":
    raise SystemExit(main())
