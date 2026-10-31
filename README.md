# AoiChangChaos_rust

**Rust で実装するカオス系データ同化。カルマンフィルタから融合粒子フィルタまで**

![Rust](https://img.shields.io/badge/Rust-2018-000000?logo=rust)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

40次元 Lorenz-96 のカオス軌道にガウスノイズを加え、Kalman Filter（KF）、Ensemble Kalman Filter（EnKF）、Merging Particle Filter（MPF）の3手法で状態推定を行うデータ同化の個人プロジェクトです。
フィルタ本体、Lorenz-96 の4次ルンゲクッタ積分、観測生成、RMSE 評価、gnuplot 用データ出力までを Rust workspace として実装しています。

> [!NOTE]
> 対象は古典的な3次元 Lorenz アトラクタではなく、40次元の **Lorenz-96** です。
> README では実装に合わせて Lorenz-96 と表記します。

## 3行で

- KF、EnKF、MPF のフィルタ本体を Rust と ndarray で実装した。
- 40次元 Lorenz-96 と観測ノイズを用意し、同一の真値に対して手法ごとの推定性能を比較した。
- 全次元観測では MPF の RMSE が 0.62、0.2 間隔の疎観測では EnKF が 2.03 となり、観測密度で手法の挙動が大きく変わった（いずれも40次元平均）。

## 技術スタック

| 領域 | 使用技術 |
| --- | --- |
| 言語 | Rust 2018 |
| 数値線形代数 | ndarray 0.14、ndarray-linalg 0.13（Intel MKL static） |
| 乱数 | rand 0.8、rand_distr、rand_xoshiro、rand_mt |
| 数値計算 | 4次ルンゲクッタ、Kalman Filter、EnKF、Merging Particle Filter |
| 検証 | 次元別 RMSE、シード固定の再現実験、出力データの比較 |
| 出力 | SSV、gnuplot |

## レビュー観点

- `docs/learning.md` に、周辺化から KF の予測と更新の式、逆行列補題、平方完成までを手で追った導出を残している。
- `predict` にクロージャを渡す設計により、EnKF と MPF が同じ Lorenz-96 の1ステップ積分を共有する。
- 出力済みの実験データと RMSE 評価器があり、結果を数値で確認できる。
- 動画シリーズとの差異を切り分けた過程を `docs/results.md` に残している。

## 対象モデル

Lorenz-96 は、循環添字を持つ次の40次元常微分方程式で表されます。

$$
\frac{dx_i}{dt} = (x_{i+1} - x_{i-2})x_{i-1} - x_i + F,
\qquad i = 0, 1, \dots, 39
$$

実験では $F = 8.0$、時間刻み $\Delta t = 0.01$ とし、4次ルンゲクッタ法で10000ステップ（$t = 0$ から $100$）積分します。
初期値は全次元を $8.0$ にし、$x_0$ だけに $10^{-3}$ の摂動を与えます。
このわずかな差が時間とともに拡大することが、データ同化を必要とする理由です。

観測は次式で生成します。

$$
y_t = x_t + w_t, \qquad w_t \sim \mathcal{N}(0, I_{40})
$$

![Lorenz-96 の軌道](pictures/1e_minus3.png)

*Lorenz-96 の軌道を3次元に投影したもの。初期値は全次元 $8.0$、$x_0$ のみ $8.0 + 10^{-3}$。*

![Lorenz-96 の観測](pictures/obs_1e_minus3.png)

*同じ軌道に標準偏差1のガウスノイズを加えた観測値。*

1次元の予備実験では、次の sin 波に標準偏差0.5のガウスノイズを加えたデータを使います。

$$
y_t = \sin(t) + w_t, \qquad w_t \sim \mathcal{N}(0, 0.5^2)
$$

## 実装したフィルタ

### Kalman Filter

`kalman_filter` は、線形ガウス系に対する厳密なベイズフィルタです。
状態遷移行列 $F$、制御入力行列 $G$、観測行列 $H$、システムノイズ共分散 $Q$、観測ノイズ共分散 $R$ を保持し、予測と更新を行います。

$$
x_{t|t-1} = F x_{t-1|t-1}, \qquad
P_{t|t-1} = F P_{t-1|t-1} F^\top + GQG^\top
$$

$$
K_t = P_{t|t-1}H^\top \left(HP_{t|t-1}H^\top + R\right)^{-1}
$$

$$
x_{t|t} = x_{t|t-1} + K_t\left(y_t - Hx_{t|t-1}\right), \qquad
P_{t|t} = P_{t|t-1} - K_tHP_{t|t-1}
$$

1次元 sin の実験では、ランダムウォークモデル $x_t = x_{t-1} + v_t$ を仮定し、$F = G = H = I$、$Q = 1$、$R = 10$ としています。

### Ensemble Kalman Filter

`ensemble_kalman_filter` は、状態分布を $N$ 個のアンサンブルで近似するフィルタです。
非線形な1ステップ予測をクロージャとして受け取り、各メンバーを独立に発展させます。

アンサンブル平均 $\bar{x}$ とサンプル共分散 $P^e$ は次式で計算します。

$$
P^e = \frac{1}{N-1}\sum_{i=1}^{N}\left(x_i - \bar{x}\right)\left(x_i - \bar{x}\right)^\top
$$

$$
K = P^e H^\top \left(HP^eH^\top + R^e\right)^{-1}
$$

$$
x_i^a = x_i^f + K\left(y + \epsilon_i - Hx_i^f\right),
\qquad \epsilon_i \sim \mathcal{N}(0, R)
$$

観測値に摂動 $\epsilon_i$ を加えることで、更新後のアンサンブル分散が過小になりにくくしています。
sin の実験では1000メンバー、Lorenz-96 の実験では300メンバーを使います。

### Merging Particle Filter

`merging_particle_filter` は、粒子フィルタの重み付きリサンプリングに「融合」を組み合わせた手法です。
各粒子の尤度を計算し、正規化した重み $\beta_i$ に従ってリサンプリングします。

$$
p(y \mid x_i) \propto
\exp\left[-\frac{1}{2}\left(y - Hx_i\right)^\top R^{-1}\left(y - Hx_i\right)\right]
$$

$$
\beta_i = \frac{p(y \mid x_i)}{\sum_{j=1}^{N} p(y \mid x_j)}
$$

リサンプリング後、$M$ 個の粒子を重み $\alpha_m$ で融合します。

$$
x_i^{\text{merge}} = \sum_{m=1}^{M}\alpha_m x_{j_m},
\qquad \sum_{m=1}^{M}\alpha_m = 1
$$

実験では $M = 3$、$\alpha = \left(\frac{3}{4}, \frac{\sqrt{13}+1}{8}, -\frac{\sqrt{13}-1}{8}\right)$ としています。
粒子数は sin と Lorenz-96 のどちらも1000です。
負の重みを含みますが、和は1になります。
融合は、リサンプリングで失われる粒子の多様性を補うための処理です。

尤度の計算では、最大対数尤度を引いてから指数関数に通すことで、アンダーフローを防いでいます。

## リポジトリ構成

```
.
├── Cargo.toml                     # workspace 定義
├── kalman_filter/                 # KF 本体と sin の example
├── ensemble_kalman_filter/        # EnKF 本体と sin / Lorenz-96 の example
├── merging_particle_filter/       # MPF 本体と sin / Lorenz-96 の example
├── lorenz96/                      # Lorenz-96 の RK4 積分とリアプノフ指数の example
├── lorenz96_generator/            # Lorenz-96 に観測ノイズを付与するバイナリ
├── sin_generator/                 # sin 波に観測ノイズを付与するバイナリ
├── rmse_evaluator/                # 真値と推定値から次元別 RMSE を計算するバイナリ
├── read_and_write/                # SSV 形式の読み書きライブラリ
├── docs/
│   ├── learning.md                # KF の導出ノート
│   └── results.md                 # 実験ログと動画との比較
├── outputs/
│   ├── Season1/                   # sin の実験結果
│   └── Season2/                   # Lorenz-96 の実験結果
└── pictures/                      # README と実験ログ用の図
```

各フィルタは独立したクレートになっています。
`read_and_write` を共通の入出力層として使い、フィルタ同士は直接依存しません。
`ensemble_kalman_filter` と `merging_particle_filter` は、Lorenz-96 の example で `lorenz96` クレートを使います。

## 実験結果

### 1次元 sin

| KF | EnKF | MPF |
| --- | --- | --- |
| ![KF](pictures/input_kf_x.png) | ![EnKF](pictures/input_enkf_x.png) | ![MPF](pictures/input_mpf_x_6.png) |

KF は観測ノイズを平滑化し、真の sin 波に近い軌道を復元しています。
EnKF と MPF も同じ観測列に対して同程度に追従しており、1次元の線形に近い問題では3手法の差が小さいことを確認できます。

### Lorenz-96

RMSE は、次元ごとの時系列 RMSE を40次元で平均した値です。

$$
\mathrm{RMSE}_d = \sqrt{\frac{1}{T}\sum_{t=1}^{T}\left(\hat{x}_{d,t} - x_{d,t}\right)^2}
$$

| 実験 | 観測する次元 | 観測間隔 | EnKF | MPF |
| --- | ---: | ---: | ---: | ---: |
| 全次元を毎ステップ観測 | 40 / 40 | 0.01 | 0.685 | 0.617 |
| 全次元を0.2間隔で観測 | 40 / 40 | 0.20 | 2.030 | 4.890 |
| 偶数次元のみ観測 | 20 / 40 | 0.20 | 3.598 | - |
| 前半20次元のみ観測 | 20 / 40 | 0.20 | 3.366 | - |

現在の example は、EnKF が全次元を毎ステップ、MPF が全次元を0.2間隔で観測する設定です。
表中の他の条件は、`observation_span`、`observation_index` の増分、$H$ を切り替えて生成しました。
数値は `outputs/Season2/part6` の RMSE ファイルから算出しています。

観測値そのものの真値に対する RMSE は、40次元平均で 1.002 です。
全次元を毎ステップ観測した場合、EnKF と MPF はともに観測ノイズ以下の RMSE を達成しています。

![MPF の全次元観測 RMSE](pictures/true_mpf_res_all_rmse.png)

*MPF で全次元を観測したときの次元別 RMSE。40次元の平均は 0.617。*

![EnKF のアンサンブル分散](pictures/enkf_mv_all.png)

*EnKF の0次元目におけるアンサンブル分散。全次元観測では分散が約5で維持され、アンサンブルが潰れていないことを確認できる。*

![EnKF の疎観測 RMSE](pictures/true_enkf_res_skip02_rmse.png)

*EnKF で0.2間隔の疎観測を行ったときの次元別 RMSE。40次元の平均は 2.030。*

![MPF の疎観測 RMSE](pictures/true_mpf_res_skip02_rmse.png)

*MPF で0.2間隔の疎観測を行ったときの次元別 RMSE。40次元の平均は 4.890。*

![EnKF の部分観測 RMSE](pictures/true_enkf_res_skip02_each_and_first_half_rmse.png)

*EnKF で偶数次元のみを観測した場合（Each）と前半20次元のみを観測した場合（First-half）の次元別 RMSE。*

全次元観測では MPF の RMSE は EnKF をわずかに下回ります。
一方、観測間隔を0.2に広げると、MPF の誤差は 4.890 まで増加し、EnKF の 2.030 を大きく上回ります。
観測が疎になると、観測時点の事前分布が広がって尤度が相対的に鋭くなり、重みが少数の粒子に集中しやすくなります。
MPF はその状態でリサンプリングするため、粒子の多様性を失いやすいことが一因と考えられます。
EnKF はアンサンブル共分散を通じて未観測方向の相関を反映できるため、疎観測でも推定が破綻しにくいという結果になりました。

前半20次元のみを観測した場合、観測した次元の RMSE は約2、未観測の次元は約4.5 になり、観測次元と未観測次元の差がはっきり出ます。
偶数次元のみを観測した場合は、観測の有無が次元ごとに交互に変わるため、RMSE が全次元で3から4程度に均されます。

## 使い方

### 必要なもの

- Rust のツールチェーン（cargo）
- Intel MKL（`ndarray-linalg` の `intel-mkl-static` feature を使うため）
- gnuplot（図を再生成する場合のみ）

`ndarray-linalg` は `intel-mkl-static` feature でビルドするため、初回ビルド時に MKL の取得とリンクが走ります。
MKL を使えない環境では、各フィルタの `Cargo.toml` の feature を `openblas-static` などに変更してください。

### ビルド

```bash
git clone https://github.com/soseigocho/AoiChangChaos_rust.git
cd AoiChangChaos_rust
cargo build --release
```

### 1次元 sin の実験

```bash
# 観測データの生成（出力先、時間間隔、シミュレーション時間、乱数シード）
cargo run --release -p sin_generator -- outputs/Season1/sin_gen.dat 0.2 200 42

# KF
cargo run --release -p kalman_filter --example sin_test -- \
  outputs/Season1/sin_gen.dat outputs/Season1/kf_da.ssv

# EnKF
cargo run --release -p ensemble_kalman_filter --example sin_test -- \
  outputs/Season1/sin_gen.dat outputs/Season1/enkf_da.ssv

# MPF
cargo run --release -p merging_particle_filter --example sin_test -- \
  outputs/Season1/sin_gen.dat outputs/Season1/mpf_da.ssv
```

### Lorenz-96 の実験

```bash
# 真値の生成（10000ステップ、t = 0.01 から 100 まで）
cargo run --release -p lorenz96 --example naive -- \
  outputs/Season2/origin_plus_1e_minus3.ssv

# 観測ノイズの付与（標準偏差1）
cargo run --release -p lorenz96_generator -- \
  outputs/Season2/origin_plus_1e_minus3.ssv outputs/Season2/obs_1e_minus3.ssv

# EnKF（全次元を毎ステップ観測）
cargo run --release -p ensemble_kalman_filter --example lorenz96_try_e -- \
  outputs/Season2/obs_1e_minus3.ssv \
  outputs/Season2/part6/enkf_mv_all.ssv \
  outputs/Season2/part6/enkf_res_all.ssv

# MPF（全次元を0.2間隔で観測）
cargo run --release -p merging_particle_filter --example lorenz96_try -- \
  outputs/Season2/obs_1e_minus3.ssv \
  outputs/Season2/part6/mpf_mv_skip02.ssv \
  outputs/Season2/part6/mpf_res_skip02.ssv
```

### RMSE の計算

```bash
cargo run --release -p rmse_evaluator -- \
  outputs/Season2/origin_plus_1e_minus3.ssv \
  outputs/Season2/part6/enkf_res_all.ssv \
  outputs/Season2/part6/true_enkf_res_all_rmse.ssv \
  estimated
```

最後の引数に `estimated` を指定すると、推定値の先頭行（$t = 0$）を読み飛ばします。
推定値と真値の両方が $t = 0$ を含まない場合は、任意の文字列を指定してください。
`rmse_evaluator` は次元別の RMSE を1行ずつ出力します。表の40次元平均は、この出力の平均です。

### 出力形式

出力はスペース区切りの SSV 形式です。
1列目が時刻、2列目以降が状態ベクトルの各次元です。

```text
0.01 8.000989796370856 7.999993662893535 ...
0.02 8.000978191325563 7.999974911880485 ...
```

sin の観測ファイルには、先頭に `time_span` と `simulation_time` のヘッダが付きます。

### パラメータの変更

現在のサンプルコードは、実験パラメータを example 内の定数として持っています。
観測間隔、観測行列 $H$、観測ノイズ共分散 $R$、アンサンブル数、粒子数は、該当する example を編集して変更します。
疎観測と部分観測の実験は、`observation_span`、`observation_index` の増分、$H$ の定義を切り替えて実行しました。
MPF の0.2間隔の example では、観測ファイルが $t = 0.01$ から始まることを踏まえ、開始インデックスを19、増分を20としています。
EnKF で同じ間隔にする場合も、開始インデックスを19、増分を20にします。
全次元を毎ステップ観測する MPF では、`observation_span` を0.01、開始インデックスを0、増分を1にします。
`ensemble_kalman_filter/examples/lorenz96_try_e.rs` には、偶数次元のみと前半20次元のみの $H$ がコメントアウトされた形で残っています。

## 実装のポイント

### モデルをクロージャで注入する

EnKF と MPF の `predict` は、次のような関数を受け取ります。

```rust
F: Fn(&Array2<f64>, &Array2<f64>) -> Array2<f64>
```

第1引数がアンサンブル状態、第2引数がシステムノイズです。
Lorenz-96 の1ステップ積分を外から渡すことで、フィルタ本体を特定のモデルに依存させていません。
sin のランダムウォークと Lorenz-96 の RK4 積分を、同じ EnKF / MPF の実装で扱えます。

### 数値的に安定な粒子フィルタ

MPF の重み計算では、対数尤度の最大値を引いてから指数関数に通します。

```rust
let max_log_likelihood = log_likelihoods
    .iter()
    .fold(0.0_f64 / 0.0_f64, |m, v| v.max(m));
let ψs = log_likelihoods
    .iter()
    .map(|log_likelihood| (log_likelihood - max_log_likelihood).exp())
    .collect::<Vec<f64>>();
```

尤度が小さい場合でも、正規化の過程で0除算やアンダーフローが起きにくくなります。

### 観測時刻の整合性チェック

各フィルタのループでは、観測データの時刻とループの時刻が一致しているかを確認し、ずれが大きい場合は panic します。

```rust
if (observation_time as f64 - loop_index as f64 * Self::ΔT).abs() > Self::ΔT / 10.0 {
    panic!("inconsistent observation datum");
}
```

観測インデックスのずれは、一見するとフィルタの性能劣化に見えるため、早い段階で落とすようにしています。

### 再現性

フィルタの乱数シードはコード内で固定しています。
再現に使うシードは次の通りです。

- sin の観測生成: CLI の第4引数（実行例では `42`）
- sin の EnKF: `808`
- sin の MPF: `801`
- Lorenz-96 の EnKF / MPF: `2039`
- Lorenz-96 の観測生成: `101`

Lorenz-96 の観測生成にはメルセンヌ・ツイスタ（`rand_mt`）を使っています。
フィルタ側のアンサンブル生成には `xoshiro256**` を使っています。

## 詰まりどころと学び

### 理論と実装の対応

`docs/learning.md` には、同時分布の周辺化、一期先予測、フィルタ更新、逆行列補題、平方完成を手で追ったノートがあります。
特に、カルマンゲインの導出で出てくる逆行列補題は、式変形を一度自分で追わないと実装の形が見えませんでした。
`ndarray-linalg` の `inv()` に渡す行列の形と、転置の位置を間違えると、コンパイルは通るのに推定だけが壊れます。

### 共分散と分散の混同

観測ノイズ $R$ はスカラーではなく共分散行列です。
1次元の sin では $1 \times 1$ 行列ですが、Lorenz-96 では $40 \times 40$ 行列になります。
MPF の尤度計算では $R$ の行列式と逆行列が必要になるため、対角成分だけを見て分散として扱うと式が合いません。
`ndarray` の `det()` と `inv()` を使い、次元を明示して実装しています。

### 粒子フィルタの退化

MPF は尤度の高い粒子にサンプルが集中するため、リサンプリングを繰り返すと粒子の多様性が失われます。
融合ステップはこの退化を緩和するための処理ですが、観測が疎になると効果が限られます。
Lorenz-96 の疎観測で MPF の RMSE が EnKF より大きくなったのは、この性質が表れた結果と考えられます。

### 動画との差異の切り分け

最初は動画シリーズと RMSE が一致しませんでした。
乱数生成器の違い（xoshiro とメルセンヌ・ツイスタ）を疑って実装を切り替えたり、観測インデックスの進め方を確認したりしました。
最終的には実験を再現でき、その過程を `docs/results.md` に残しています。
結果が合わないときに、乱数、観測スケジュール、RMSE の比較方法を順に切り分けた経験は、数値実験のデバッグに役立ちました。

## 既知の制約

- 実験パラメータが example にハードコードされており、観測間隔、$H$、$R$、アンサンブル数、粒子数を CLI や設定ファイルから指定できない。疎観測と部分観測の再現にも example の編集が必要。
- フィルタの乱数シードは固定だが、単一シードの結果しかなく、複数シードの統計評価を行っていない。
- 観測生成のノイズ分散と、フィルタ側の観測ノイズ共分散が一致していない条件がある。Lorenz-96 の観測は分散1で生成しているが、EnKF は標準偏差3の摂動観測から $R$ を見積もり、MPF は $R = 50I$ としている。手法間の厳密な性能比較には $R$ の統一が必要。
- 自動テストと CI が未整備で、検証は出力済みデータと RMSE の目視確認に依存している。
- 可視化は gnuplot のスクリプトに依存しており、Rust のプログラムから直接図を出していない。
- 対象は Lorenz-96 のみで、古典的な3次元 Lorenz アトラクタは含んでいない。

## 改善案

- 観測スケジュール、$H$、$R$、アンサンブル数、粒子数を設定ファイルまたは CLI から指定できるようにする。
- 予測と更新の形状チェック、共分散の対称性、PSD 性を確認するテストを追加する。
- 複数シードで RMSE の平均と分散を評価し、手法比較を統計的に行う。
- EnKF に localization と inflation を追加し、疎観測での性能を検証する。
- MPF のリサンプリングと融合のパラメータを感度分析する。
- Rust のプロットライブラリかノートブックで可視化を完結させる。
- GitHub Actions でビルドとテストを自動化する。

## ドキュメント

- [`docs/learning.md`](docs/learning.md): カルマンフィルタの導出ノート。周辺化、一期先予測、フィルタ更新、逆行列補題、平方完成の証明を含む。
- [`docs/results.md`](docs/results.md): 実験ログ。動画シリーズとの比較、出力ファイルの説明、再現時の注意点を含む。
- [`outputs/Season1`](outputs/Season1): sin の観測データと KF / EnKF / MPF の推定結果。
- [`outputs/Season2`](outputs/Season2): Lorenz-96 の真値、観測値、リアプノフ指数、RMSE。
- [`outputs/Season2/part6`](outputs/Season2/part6): EnKF と MPF の推定結果、アンサンブル分散、RMSE。

## クレジット

本リポジトリは、のーふぉん氏（Twitter: [@4th_No_Fon](https://twitter.com/4th_No_Fon)）制作の動画シリーズ「[琴葉葵の混沌解析](https://www.nicovideo.jp/series/236665)」で紹介されたデータ同化シミュレーションを、Rust で実装し直したものです。
動画シリーズの解説と問題設定を参考にしていますが、本リポジトリは個人による独立した実装であり、動画の配布物そのものではありません。

LICENSE は [MIT ライセンス](LICENSE) です。
Copyright (c) 2021 粗製伍長。

## English Summary

A personal Rust workspace that implements data assimilation for a 40-dimensional Lorenz-96 system.
It contains from-scratch implementations of the Kalman Filter, Ensemble Kalman Filter, and Merging Particle Filter, along with RK4 integration, observation generation, RMSE evaluation, and reproducible experiment outputs.
The experiments compare the filters under dense, sparse, and partial observations, and the repository includes derivation notes and a detailed experiment log in Japanese.
