# M5: Do deep order-book models beat order-flow imbalance after costs?

_itchbook walkthrough, milestone 5. The second research study: a deep convolutional network that reads the raw order book, against a linear model on order-flow imbalance, on seven NASDAQ days with strictly out-of-time testing, scored first on forecasting and then on trades that pay the spread and NASDAQ's fees._

## The question

Short-horizon price forecasting from the limit order book is one of the most published topics in quantitative finance. DeepLOB (Zhang, Zohren & Roberts, 2019) is the best known of the deep models. It reads the top ten levels of the book over the last hundred updates and classifies the next move as down, flat or up, and it reports accuracies above 60%. The linear alternative is older and simpler. Cont, Kukanov & Stoikov (2014) showed that order-flow imbalance (OFI), the net pressure of arrivals and departures at the best quotes, explains most of the contemporaneous mid-price change, with R² around 65%. Xu, Gould & Howison (2019) extended it to ten levels.

Two questions are usually left open. The first is whether the deep model beats a well-specified linear OFI model at all, rather than a weak baseline. The second is whether either forecast survives the cost of trading on it. The literature mostly reports classification accuracy and, at best, frictionless profit at the mid-price. A forecast that is right 55% of the time is worth nothing if the move it predicts is smaller than the spread you pay to act on it. This milestone answers both questions on the same data, with every choice made before the test days were looked at.

## Data, universe and the split

The data is seven full NASDAQ TotalView-ITCH days that NASDAQ still publishes: 30 January, 27 March, 30 July, 30 August, 30 October and 30 December 2019, and 30 January 2020. They span a year, which matters. A model that has only learned one day's conditions should be caught out by days months later.

The universe was fixed by a rule written down before any data was examined. The rule takes the twenty symbols with the most trade prints on the first day, 30 January 2019, among those trading at $5 or more. The result is SPY, AMD, AAPL, MSFT, NVDA, FB, QQQ, BABA, MU, TSLA, AMZN, T, QCOM, INTC, IWM, TQQQ, CRON, EBAY, VZ and SBUX. All twenty are present on all seven days. The price floor keeps out sub-$5 names, whose book structure is quite different and which would dominate a pooled model. NASDAQ's feed includes its trading in NYSE-listed ETFs such as SPY and IWM, which is why they appear.

The split is by date and is never shuffled:

- **Training:** the first four days.
- **Validation:** 30 October 2019. Every choice is made on it: the ridge penalty, the classification cut-offs, the trading thresholds and early stopping.
- **Test:** 30 December 2019 and 30 January 2020. These are scored once, at the end.

A random split of rows would leak information, because neighbouring rows overlap in time and share the same market conditions. A date split doesn't.

## From the book to features

The L3 engine from M2 and M3 already holds the true book after every message. `itch-lob` samples it on a fixed 250-millisecond grid over regular hours, 93,600 rows per symbol per day. Each row is the state after every message up to its time stamp, the same convention as the M2 snapshots, and it records three things:

- **The book itself:** the top ten levels on each side.
- **Multi-level OFI** for levels one to ten, accumulated over every book event in the interval rather than computed between snapshots. A 250 ms snapshot misses most of what happens inside it, and OFI is about exactly those events.
- **Traded volume by aggressor side.** An execution against a resting buy order is seller-initiated.

The OFI for one event at level m follows Cont, Kukanov & Stoikov. The bid side adds the new size if the bid price rose or held, and subtracts the old size if it fell or held. The ask side does the mirror image with the opposite sign. So a growing bid queue or a shrinking ask queue counts as buying pressure. A missing level is treated as price zero on the bid side and infinity on the ask side, so a level that appears or vanishes counts as a full arrival or departure.

This code is held to the project's standard. A property test replays random event streams into the recorder and into an independent replay built on the M2 `ReferenceBook`. The replay recomputes every row from scratch: the top levels at each grid time, the OFI from the book before and after each event (with the formula written out case by case, independently), and the signed volume from each resting order's side. Three bugs were planted to check that the test notices them:

- moving the grid boundary by one nanosecond;
- flipping one sign in the ask term;
- swapping buyer- and seller-initiated volume.

Each one was caught. On real data the contemporaneous correlation between best-level OFI and the mid-price change in the same interval is 0.67 to 0.75 across symbols, in line with the published R² of about 65%. One detail of the feed showed up along the way: NASDAQ reports every hidden execution with side "B", so the side of hidden trades is not known and they are excluded from the signed volume.

The models see three views of each row.

- **The book view** is DeepLOB's input. For each of the ten levels it holds the ask price, ask size, bid price and bid size. Prices are in ticks from the mid and sizes are log round lots, so a $300 stock and a $30 stock look alike and one model can serve all twenty.
- **The order-flow view** is the ten OFI levels plus signed volume, scaled by the symbol's typical best-level depth, and the event count.
- **The linear features** are OFI and signed volume summed over the last 1, 4, 20 and 120 intervals (0.25 s to 30 s), plus book imbalance at the best level and over five levels, the spread and the event count. That is 48 features.

All scales are computed from the training days only.

## Targets and how forecasts are scored

The targets are the mid-price change over 1, 5 and 30 seconds (4, 20 and 120 grid steps), in ticks, and its sign: down, flat or up. The class balance changes sharply with the horizon:

- **1 second:** the mid doesn't move in 63% of test rows.
- **5 seconds:** 37% of rows are flat.
- **30 seconds:** 18% of rows are flat.

That imbalance is why accuracy alone misleads: always predicting "flat" scores 63% at one second. Classification is therefore scored with macro-F1, the average F1 of the three classes, where that constant guess scores 0.26.

A second number isolates direction. When the model calls a move and the price really moves, how often is the direction right?

## The models

**Ridge on multi-level OFI (the baseline).** A linear regression of the tick change on the 48 linear features, with a ridge penalty chosen on the validation day. It is fitted in a single streaming pass over the training days, accumulating XᵀX, so the training set never has to be in memory. Its validation R² is small: 0.0033, 0.0031 and 0.0008 for the three horizons. Predicting a price change from what happened before it is far harder than explaining the change from what happened during it.

To turn the forecast into a class, a forecast below a cut-off in absolute value counts as "flat". The cut-off is chosen per symbol and horizon on the validation day. A per-symbol cut-off matters for fairness: on a one-tick-spread stock like AMD the mid rarely moves, while on AMZN it almost always does. A deep model can learn that from the book, and a single shared cut-off would handicap the baseline.

**DeepLOB.** The network follows the paper's layout:

- Three convolution blocks merge each level's price with its size, then the bid side with the ask side, then all ten levels into one.
- Each block also carries convolutions over time.
- An inception module combines several time scales, an LSTM with 64 units summarises the hundred-step (25 s) window, and a three-class output head serves each horizon.

It follows the official implementation's choices: 32 filters, batch normalisation, and 64 channels per inception branch. There are two deliberate differences:

- **One network, three outputs.** It is trained jointly on all three horizons rather than as three separate networks, which saves two-thirds of the compute.
- **Faster optimiser settings.** Adam with learning rate 10⁻³ and batches of 256, where the paper uses 10⁻² with ε = 1 and batches of 32, chosen for throughput on the available hardware.

**The order-flow deep model.** The same network, with the two level-merging convolution blocks replaced by a single convolution across the twelve order-flow inputs: OFI at levels one to ten, signed volume and the event count. Kolm, Turiel & Westray (2023) and Lucchese, Pakkanen & Veraart (2024) both found that deep models do better on order-flow inputs than on raw book levels. This variant tests whether that holds here. With twelve inputs instead of forty it trains far faster, up to about 2,000 windows a second on the same laptop.

## From forecast to trade

A forecast earns money only if acting on it costs less than it predicts. The trading test is deliberately simple and identical for every model.

- **When it decides.** At decision points one horizon apart, so positions never overlap, the model's expected move is compared with the round-trip cost of a taker trade.
- **The cost.** Cross the spread now (half a spread from the mid), cross it again on exit (assumed to be the same spread), and pay NASDAQ's fee for removing liquidity twice. That fee was $0.0030 per share in 2019–20, per the exchange's price list in its SEC filings.
- **The decision.** If the expected move beats the cost by more than a threshold, the model buys at the ask or sells at the bid, and unwinds at the opposite touch one horizon later.
- **Expected moves.** For ridge, the expected move is its forecast in ticks. For the deep models, it is each class's probability times that symbol's average move in that class on the training days.
- **The threshold.** Chosen on the validation day to maximise total profit after costs. Not trading, which earns exactly zero, is always one of the options, so a model with no usable edge simply doesn't trade on the test days.

Two references bracket the results:

- **Frictionless.** The same decisions settled at the mid with no spread and no fee. It measures the raw forecasting edge in ticks.
- **The oracle.** It knows the future mid-price and trades only when the realised move beats the cost. It is the upper bound on what any forecast could earn under these costs.

## Results

On the two test days, scored once:

| | Ridge on OFI | DeepLOB | Deep model on OFI |
|---|---|---|---|
| Macro-F1, 1 s | 0.426 | **0.465** | 0.435 |
| Macro-F1, 5 s | 0.451 | **0.533** | 0.526 |
| Macro-F1, 30 s | 0.426 | **0.505** | 0.490 |
| Accuracy, 1 s (majority class: 0.633) | 0.564 | **0.663** | 0.660 |
| Direction right when it calls a move, 1 s | **62.7%** | 55.7% | 53.6% |
| Direction right when it calls a move, 30 s | 52.4% | 52.1% | 51.6% |

DeepLOB is the better classifier at every horizon, beating the linear model by 0.04 to 0.08 in macro-F1. That is consistent with the literature. Where the gain comes from is the important part. DeepLOB's accuracy at one second is ten points higher, but when it calls a move, it gets the direction right *less* often than the linear model: 55.7% against 62.7%. Most of its advantage is in recognising when nothing will happen. That is valuable for classification scores and of no use to a trader, who only acts on calls of direction.

The order-flow deep model sits between the other two. At five seconds it nearly matches DeepLOB (0.526 against 0.533). At one second it is barely above the linear model (0.435 against 0.426), and its direction calls are the least accurate of the three. So here, unlike in Kolm et al. and Lucchese et al., order-flow inputs did not beat book inputs on classification. They did give the larger edge before costs at one and five seconds (see below).

The trading test makes that concrete:

![Edge per trade before costs, mid to mid, on the test days, against the cost of a round trip on a stock with a one-tick spread with and without the taker fee.](docs/figures/m5_edge.png)

| Ticks per trade, test days | 1 s | 5 s | 30 s |
|---|---|---|---|
| Ridge, frictionless | 0.14 | 0.25 | 0.23 |
| DeepLOB, frictionless | 0.11 | 0.21 | 0.43 |
| Deep model on OFI, frictionless | 0.14 | 0.33 | 0.39 |
| Any model, after the spread (with or without the fee) | no profitable trades | no profitable trades | no profitable trades |
| Oracle, after the spread and the fee | 1.1 | 2.6 | 7.6 |

All three models have a real edge before costs: on average the mid-price moves their way, by about a tenth to four-tenths of a tick. No model leads everywhere. The order-flow deep model has the largest edge at five seconds (0.33 ticks) and DeepLOB at thirty (0.43). Most of these stocks trade at a one-tick spread most of the time. The exceptions on the test days are AAPL, FB, BABA, NVDA, TSLA and AMZN, plus MSFT, QCOM and TQQQ on one of the two days. So the cheapest possible round trip costs one tick of spread, and 1.6 ticks with the fee. The edge is several times smaller than the cheapest cost. With the thresholds chosen on the validation day and the spread alone as the cost, the rule finds no trades worth taking at one and five seconds. At thirty seconds the handful it takes lose money: the ridge rule took 141 trades and lost 2.3 ticks per trade, DeepLOB took 3 and the order-flow model none. With the fee added, the ridge rule took 9 and 20 trades at one and five seconds and lost on both.

The oracle shows that this isn't because there is nothing to earn. Up to 7.6 ticks per trade survive the costs for anyone who knows the future. None of the models captures enough of the size of the move to clear the spread.

## Why a better forecast still can't trade

A probability of 55% on the direction is a large statistical edge, and with millions of rows it is far from chance. But the payoff of a taker trade is the move minus the spread. At these horizons the typical move is less than a tick, and the spread is at least one. A forecast must therefore predict not just the direction but large moves, and predict them reliably enough that the average winner exceeds a tick plus fees. Short-horizon book data mostly predicts the next small step, often the next change at the touch, which is exactly the move the spread swallows. DeepLOB's extra skill is concentrated in the "no move" class, which earns nothing. The thirty-second horizon is where DeepLOB's edge per trade is largest. It is also where predictability is weakest: by thirty seconds every model gets the direction right only about 52% of the time.

This agrees with what has been reported where people looked:

- **LOBCAST** (Prata and co-authors, 2024): benchmarked fifteen deep models and found their performance drops sharply on new data. Its trading backtest, run without fees, was profitable only for some stocks.
- **Briola, Bartolucci & Aste (2025):** found that forecasting power does not necessarily translate into actionable signals.

What this study adds is the same verdict reached with a fair linear baseline, a strictly out-of-time split, and costs from the exchange's own price list.

## Limits of the study

- **Training on a laptop.** Both deep models were trained on an 8 GB laptop, in passes of 256,000 windows sampled from the 7.5 million training windows. A first run was cut short by memory pressure. A second run gave each model more patience before stopping: training stopped only after five passes in a row without improvement on the validation day. DeepLOB stopped after 11 passes (2.8 million windows). Its validation score had peaked at the sixth pass, the same model the first run had kept, so its results did not change; training is seeded and exactly reproducible. The order-flow model also stopped after 11 passes. Its best pass improved its test scores slightly over the first run (macro-F1 up by 0.003 to 0.023, edge before costs up to 0.39 ticks), and those are the numbers reported here. Validation performance had levelled off for both. Bigger models or many more days of data might still help. No model's edge comes within a factor of three of the cost of a round trip, so neither change is likely to alter the trading verdict.
- **Two test days** are a small test set for the trading result, and they are not a market-regime study.
- **Only taking liquidity is tested.** A passive strategy that posts limit orders earns the spread instead of paying it, but then the question becomes whether it is filled. M4 showed that L2 models of queue position, which such a backtest would need, are biased in exactly the cases that matter. Testing passive execution properly needs the L3 queue simulation this project's engine provides, which is a natural next study.
- **The grid is 250 ms in clock time,** not event time. Some very fast signal is averaged away.
- **Results pool twenty large, liquid names.** Small stocks with wide spreads behave differently.

## Questions to expect

*Why compare against ridge on OFI instead of a naive baseline?* Because a deep model that beats a constant guess proves very little. OFI is the established linear description of short-term price pressure, and a ten-level version with several look-back windows is a strong, cheap competitor. If a deep model can't clearly beat it, its extra complexity isn't earning anything.

*Why does DeepLOB win on macro-F1 but lose on direction?* Macro-F1 rewards getting the "flat" class right, and DeepLOB is much better at predicting when the mid won't move. On the rows where it commits to a direction and the price does move, the linear model is right more often. A trader only acts on those rows.

*Isn't "no profitable trades" just a bad threshold?* The threshold is chosen on a separate day to maximise profit after costs, and "don't trade" is always allowed. The frictionless results show the edge directly: 0.1 to 0.4 ticks per trade. That is below the cost of one round trip at a one-tick spread, before any fee.

*What would change the verdict?* Execution without paying the spread, meaning passive orders, which need a correct fill model (see M4). Or a much longer horizon with more predictable moves. Or a cost structure with large maker rebates. The study's code supports all three as extensions.
