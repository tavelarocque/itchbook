# M4: How wrong are L2 queue-position models?

_itchbook walkthrough, milestone 4. The first research study: what a market-by-price observer can know about its place in line, how the standard estimation models behave against the true L3 queue on two NASDAQ days, and where the error comes from._

## The question

A resting limit order earns the spread only if it is filled, and on a FIFO exchange it is filled only when everything ahead of it at its price has gone. So the single most important number for a passive strategy is the number of shares ahead of its order. A strategy that sees the full order-by-order (L3) feed can count them. Most researchers and most backtests don't have L3 data. They have L2: the total size at each price, plus the trades. From L2 they must estimate the queue ahead, and every backtest of a passive strategy rests on that estimate. If the estimate is too optimistic, the backtest fills orders that would never have been filled, and the strategy looks better than it is. If it is too pessimistic, the opposite happens.

The difficulty is a single ambiguity. When the size at a price falls and there is a trade, the trade came off the front of the queue, so the order moved up. When the size falls and there is no trade, someone cancelled, and L2 doesn't say whether that someone was ahead of the order or behind it. Every L2 queue model is an answer to that one question. With the L3 book from M2 and M3, the answer can be measured instead of assumed. This milestone measures it.

## What the observer sees, and the models

For each order being followed, the study gives the models exactly what an L2 feed would show about the rest of its price level: the total size of the other orders there, message by message, and every trade print at that price. The order's own size is excluded, because a backtest that simulates an order never sees it in the data. The models are not told when the order itself is filled. They see those prints like any other.

A trade of v shares always comes off the front, so every model reduces its estimate of the shares ahead, A, by v (not below zero). If the print is larger than A, the model concludes the order was reached and predicts a fill. For a cancellation of c shares the models differ. Write B for the shares behind. Each model says what fraction φ of the cancelled shares came from ahead:

    A' = A − c · φ(A, B)

The result is clamped to what is possible: at most c shares can leave from ahead, and at most B from behind. The models tested are these.

- **Pessimistic** (φ = 0): cancels always come from behind, so only trades move the order up. This is hftbacktest's `RiskAdverseQueueModel`, the conservative default in practice.
- **Optimistic** (φ = 1): cancels come from ahead first.
- **Uniform** (φ = A/(A+B)): every share in the queue is equally likely to be cancelled. It is the pro-rata assumption behind most simple estimators and several academic models.
- **power2, power3, log**: hftbacktest's `ProbQueueModel` with `PowerProbQueueFunc(2)`, `PowerProbQueueFunc(3)` and `LogProbQueueFunc`. A unit test transcribes hftbacktest's update rule and checks these give the same numbers. One deliberate difference: here an estimate never goes below zero on a cancellation.
- **front3** (φ = r³, with r = A/(A+B)): hftbacktest's `PowerProbQueueFunc3(3)`. hftbacktest's own Level-3 example compares against this one. It sends most cancellations behind.
- **Calibrated**: φ is the measured share of cancelled volume that was ahead, as a function of r. It is fitted on the *other* trading day, so it is always scored out of sample.

## Ground truth, and why it can be trusted

The true number ahead is kept incrementally, without walking queues. Every order carries its arrival sequence number. Orders only ever join a queue at the back, so an order X is ahead of a tracked order T exactly when they share a price and X arrived first. The truth for T therefore changes only when an earlier order at T's level loses shares, and it changes by exactly those shares. Replaces need care. A replace at the same price is a new order and goes to the back: in truth the old order leaves from wherever it was, even though L2 sees only the net change in size.

This is checked three ways, following the project's rule that every layer is checked against something it can't fake. A property test feeds random event streams to the study and to `FastBook`. After every event it compares each tracked order's truth with the shares ahead `FastBook` finds by walking the actual linked queue. Three bugs were planted to check that the test notices them, and it caught each one:

- counting large cancellations from behind as if they were ahead;
- letting a same-price replace keep its priority;
- counting executions behind the order as ahead of it.

An end-to-end test runs the real tool on a synthetic day with verification at every evaluation point. Finally, on both real days the tool checked its truth against the queue walk at every single evaluation point: 596,688 checks on 30 December 2019 and 801,382 on 30 January 2020, with no failures and no orphan events.

## Who was followed

The population is orders that join the best bid or ask during regular hours, while the stock is trading, behind at least one other order. Orders that set a new best price have nothing ahead, and orders deeper in the book are not yet in contention. A deterministic hash of the order reference picks 1 in 16 of them. That gave 4,608,261 tracked orders in 7,136 symbols on 30 January 2020, and 2,876,440 in 7,062 symbols on 30 December 2019. A whole day runs in about 95 seconds on the laptop.

The population is dominated by fleeting orders:

- The median tracked order lived 0.15 seconds.
- 90% ended by being deleted, 6% by a replace, and only 4.5% were ever filled.
- Only 6% saw even one trade at their price.

The median order joined behind 4 others and 700 shares. When an order was filled, it was at the head of its queue 97.5% of the time, as strict time priority predicts.

Because most orders die before anything happens, medians over all orders are uninformative: a zero error mostly means nothing has changed yet. The estimation error below is therefore averaged over each order's evaluation points, which are the trade prints at its price up to its first real fill. Each order's error is divided by the queue it joined behind, so a 10-share error behind 100 shares and a 1,000-share error behind 10,000 count the same. The fill results use every order, because predicting a fill for an order that was actually cancelled is exactly the error a backtest makes.

## Finding 1: cancellations come from behind

For every cancellation at a tracked order's level, the study records the order's true relative position r = ahead/(ahead + behind) and whether the cancelled shares were ahead of it or behind. The share of cancelled volume that came from ahead, as a function of r, is the true φ:

![The measured share of cancelled volume that was ahead of a tracked order, by its true position in the queue, on both days, against the uniform model, hftbacktest's power2 and its front3 (r³). The two days lie on top of each other.](docs/figures/m4_phi.png)

The curve is nowhere near uniform. An order with half the queue ahead of it (r = 0.5) would, under the uniform model, expect half of every cancellation to come from ahead. In fact, for orders between r = 0.45 and r = 0.55, 4 to 5% did, on both days. Even at r = 0.9, with nine shares ahead for every one behind, only about half of the cancelled volume came from ahead. Across all positions, 15% of cancelled volume came from ahead of the tracked order and 85% from behind. The two days differ by at most about two percentage points anywhere along the curve, and by less than one over most of it. That stability is what makes it legitimate to fit the calibrated model on one day and test it on the other.

The reason is the population of orders itself. In a FIFO queue, position and age are the same thing: the orders behind you are younger than you. Young orders are overwhelmingly the ones that get cancelled, most of them within a fraction of a second of arriving, as the 0.15-second median life shows. This agrees with what has been measured elsewhere:

- Gu and co-authors (2013) found cancellations almost absent at the head of the queue on the Shenzhen exchange.
- Dahlström, Hagströmer and Nordén (2024) found that volume behind an order, not ahead of it, drives its cancellation hazard.
- Databento's 2025 note on reconstructing queue position states that cancels are mostly pulled from the back.

None of these measures the curve against the shares ahead on US ITCH data, which is what the figure adds. Even hftbacktest's most conservative probabilistic function, r³, puts several times too much cancellation ahead in the middle of the queue.

## Finding 2: part of the queue ahead moves invisibly

The study also breaks down where the shares ahead of tracked orders actually went. Over 30 January, summed across tracked orders:

| Cause | Shares that left from ahead | From behind |
|---|---|---|
| trades | 114 M | 8 M |
| partial cancels | 148 M | 811 M |
| deletes | 3,589 M | 22,730 M |
| replaced to another price | 838 M | 2,688 M |
| replaced at the same price | 337 M | 788 M |

Same-price replaces moved 337 million shares from ahead of tracked orders to the back of their queues. That is three times more than trades did. An L2 feed sees a same-price replace only as the net change in the level's size, and for those replaces the net change totalled just 4.9 million shares. So about 6.6% of everything that left from ahead left invisibly. To an L2 observer the level size barely moves, and no model built on L2 can see it, however well calibrated. This puts a floor under the error of every L2 model. It explains why even the pessimistic and calibrated models, which get the cancellation direction nearly right, still think shares are ahead at the moment a real fill arrives.

## Finding 3: estimation error grows with the queue

The headline numbers on 30 January 2020 (30 December 2019 is within about one percentage point on every figure):

| Model | Mean error / queue at entry | Bias | Orders within 10% |
|---|---|---|---|
| pessimistic | 0.076 | +0.070 | 82% |
| optimistic | 0.146 | −0.146 | 66% |
| uniform | 0.118 | −0.092 | 67% |
| power2 | 0.105 | −0.057 | 68% |
| power3 | 0.100 | −0.041 | 69% |
| log | 0.123 | −0.102 | 69% |
| front3 | 0.081 | +0.017 | 77% |
| calibrated (fitted on 30 Dec) | 0.075 | +0.038 | 81% |

A positive bias means the model thinks more shares are ahead than really are. The uniform model, the textbook assumption, is wrong by 12% of the queue on average and always in the optimistic direction. The only models that do better are the ones that send cancellations behind: pessimistic, front3 and calibrated. Weighting every symbol equally gives the same ranking, and the top ten symbols hold under 7% of the evaluated orders, so no single stock drives the result. The ranking also holds for orders that rested at least one second.

The error depends strongly on how long the queue was:

![Median relative error of the shares-ahead estimate by the number of orders ahead when the order joined, 30 January 2020.](docs/figures/m4_error_by_queue.png)

Behind one to four orders, every model is almost always right, because there is little to get wrong. Behind 25 to 99 orders, the uniform model's mean error is 18% of the queue and the optimistic model's 25%, while the calibrated model stays at 6%. Following orders as trades arrive shows the bias directly. After every print, the uniform family sits well below the truth (black) and the pessimistic model sits above it:

![Mean shares ahead, as a share of the queue at entry, after k trade prints at the order's price, for the sampled orders' evaluation points on 30 January 2020. The rise after about ten prints is survivorship: orders that stay that long are the ones with large queues still ahead.](docs/figures/m4_trajectory.png)

## Finding 4: what this does to a backtest

A backtest needs fills, not estimates. For each tracked order the study records whether each model would have filled it while it was alive, and when. Two numbers matter:

- **Precision** is the share of predicted fills that really happened. A predicted fill that never happened is a phantom fill, profit the strategy never earns.
- **Recall** is the share of real fills the model predicted.

Aggregate fill rates hide everything. The uniform model predicts a fill rate of 4.6% against a true 4.5% and looks almost perfect, but only because its phantom fills and its missed fills roughly cancel. Precision and recall separate them:

![Precision against recall of each model's simulated fills, all tracked orders (left) and orders that joined behind 10 or more orders (right), 30 January 2020.](docs/figures/m4_fill.png)

For orders that joined behind 10 to 24 others, which is where a passive strategy's fills are actually decided:

| Model | Precision | Recall | Filled shares, predicted / real | Predicted early |
|---|---|---|---|---|
| pessimistic | 98.9% | 43.7% | 0.58 | 0.5% |
| optimistic | 71.6% | 92.2% | 1.29 | 30.5% |
| uniform | 83.4% | 85.7% | 1.10 | 24.4% |
| power2 | 91.9% | 80.8% | 0.96 | 18.4% |
| front3 | 96.8% | 63.1% | 0.75 | 8.6% |
| calibrated | 97.5% | 57.1% | 0.69 | 5.6% |

There is a clear trade-off and no model escapes it. Optimistic-leaning models catch most real fills, but one in six uniform fills and more than one in four optimistic fills never happened, and a quarter to a third of their fills arrive before the real one. A backtest using them fills orders in exactly the situations where real orders are still waiting, which is when the price is about to move through them. The pessimistic model is almost never wrong when it says "filled", but it misses more than half the real fills, because it never credits the invisible advances of Finding 2 or the few cancellations that do come from ahead. The calibrated model, fitted on the other day, has the best estimation error and 97.5% precision, but still misses 43% of fills in long queues, for the same reason.

## What this means

For a backtest on L2 data, the evidence argues for a conservative model and for reporting its error. The uniform and log models are optimistic in a systematic way that grows with queue length, and that is the direction that flatters a strategy. Fitting φ to L3 data from the same venue on a different day removes most of the estimation bias and transfers between days. But no L2 model recovers the fills lost to same-price replaces, so a conservative L2 backtest will under-fill by a predictable amount: in queues of 10 to 24 orders, the front3, calibrated and pessimistic models fill 25%, 31% and 42% fewer shares than were really filled. Where L3 data exists, it should be used directly. That is the case for this project's engine, and for M5.

## Limits of the study

- **Two trading days.** One of them is the quiet session between Christmas and New Year. The agreement between them is reassuring, but it is not a market regime study.
- **The population is orders joining the best price.** The φ curve is measured at the levels of tracked orders and is share-weighted. A level with several tracked orders counts its cancellations once per tracked order.
- **Errors are measured at trade prints.** Those are the moments a fill is decided, but the error between prints isn't measured.
- **The calibrated model applies φ at its own estimated position,** not the true one, so it inherits some of its own error.
- **The order is in the data.** Each order is followed in a world where it really was in the queue. A simulated order that wasn't there would not change the prints others received, except for the shares that would have traded against it. That is the standard counterfactual limit of any backtest, and it is the same for every model.

## Questions to expect

*Why measure against L3 rather than just assume a model?* Because the models disagree by a factor of two in the cases that matter, and they cannot all be right. On NASDAQ the L3 feed makes the truth observable for every order, so the comparison needs no assumptions about the market.

*Why is cancellation so concentrated at the back?* Position in a FIFO queue is age, and young orders are the ones that are cancelled. Most are placed and withdrawn within a fraction of a second, which the 0.15-second median life shows directly. An order that has survived long enough to be near the front belongs to a different population: patient liquidity that is rarely pulled.

*If cancellations are mostly behind, why isn't the pessimistic model right?* For two reasons. Some cancellations do come from ahead, 15% of the volume overall. And same-price replaces move orders from ahead of you to the back without changing the level size an L2 feed shows. On 30 January that was 337 million shares from ahead, three times what trades removed.

*How do you know the ground truth is right?* It is checked against the L3 book's own queue walk at every evaluation point on both full days, about 1.4 million checks with no failures. A property test does the same after every event of random streams, and that test has been shown to catch three deliberately planted bugs.

*Isn't the calibrated model just fitting the answer?* It is fitted on one day and scored on the other, in both directions, and the φ curves of the two days agree to within a percentage point.
