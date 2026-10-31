## 周辺化
ある確率変数に対して独立な確率変数との同時確率を、その独立な確率変数について積分を行った結果はその独立な確率変数を用いない場合の確率に等しい。また同時確率にベイズの定理を用いることもできる。(ただし独立性は周辺化自体には要請されない。)
$$
p(a) = \int p(a,b)db = \int p(a|b)p(b)db = \int p(b|a)p(a)db
$$
多変数ベイズは変数が独立であるならば結局事前条件が増えただけで、事前条件から引っ張り出した変数と元々の確率変数による残された事前条件における同時確率/事前条件から引っ張り出した変数による残された事前条件における確率になる。
または
確率変数を事前条件に押し込んだ確率
*押し込んだ確率変数による押し込む前の事前条件における確率
となる

補正は観測yによる推定状態xの補正(状態xは観測できないため)

where
$$
\begin{align}
    \text{推測: }x_t &= F_t x_{t-1} + G_t v_t \\
    \text{観測: }y_t &= H_t x_t + w_t
\end{align}
$$

$p(x_{t-1}|y_{1:t-1})$ -(一期先予測)-> $p(x_t|y_{1:t-1})$ -(フィルタ)-> $p(x_t|y_{1:t})$

## 一期先予測
$$
\begin{align}
    p(x_t|y_{1:t-1}) &= \int p(x_t,x_{t-1}|y_{1:t-1}) dx_{t-1} \\
        &= \int p(x_t|x_{t-1},y_{1:t-1}) p(x_{t-1}|y_{1:t-1}) dx_{t-1} \\
        &= \int p(x_t|x_{t-1}) p(x_{t-1}|y_{1:t-1}) dx_{t-1} \\
    &\because
        \substack{
            \text{
                シミュレーション$f_t$の引数に観測$y$は用いられていない
            } \\
            \text{
                (し、シミュレーションの結果に対して以前の観測は独立である)ため
            }
            \text{
                同時確率から取り除くことができる
            }
        }
\end{align}
$$

## フィルタ
$$
\begin{align}
    p(x_t|y_{1:t}) &= p(x_t|y_{1:t-1},y_t) \quad \text{スライスの分解} \\
        &= \frac{p(x_t,y_t|y_{1:t-1})}{p(y_t|y_{1:t-1})} \quad \text{ベイズ} \\
        &= \frac{p(y_t|x_t,y_{1:t-1}) p(x_t|y_{1:t-1})}{p(y_t|y_{1:t-1})} \quad \text{分子にベイズ} \\
        &= \frac{p(y_t|x_t) p(x_t|y_{1:t-1})}{p(y_t|y_{1:t-1})} \\
    &\because
        \substack{
            \text{$y_t$はある状態$x_t$について過去の観測$y_{1:t-1}$に依存しないことから} \\
            \text{事前条件から削除} \\
        } \\
        &= \frac{p(y_t|x_t) p(x_t|y_{1:t-1})}{\int p(x_t,y_t|y_{1:t-1}) dx_t} \quad \text{これは}x_t\text{について周辺化} \\
        &= \frac{p(y_t|x_t) p(x_t|y_{1:t-1})}{\int p(y_t|x_t,y_{1:t-1})p(x_t|y_{1:t-1}) dx_t} \quad \text{分母にベイズ}
\end{align}
$$

## 一期先予測の導出
推測が完璧に正しいと仮定 -> $x_t$について厳密に等号が成立 -> 次のデルタ関数に関する等号が成立
$$
p(x_t|x_{t-1},v_t) = δ(x_t - F_t x_{t-1} - G_t v_t)
$$
$\because$δ関数の引数に当たる式が0となるある$x_t$の時に確率が1に、それ以外の場合は0になるため
$$
    \sim \frac{1}{(2\pi)^{\frac{k}{2}}} \frac{1}{\sqrt{|\varepsilon^{2}I_K|}} \exp{(-\frac{1}{2} (x_t-F_t x_{t-1} - G_t v_t)^{\text{transpose}}
        \times (\varepsilon^{2}I_k)^{-1} (x_t - F_t x_{t-1} - G_t v_t))}
$$

$\because$とはいえシミュレーションは厳密でないのでδ関数を正規分布で近似

$$
\begin{align}
    &= \mathcal{N}(G_t v_t, \varepsilon^{2I_k})
\end{align}
$$

$\because$この正規分布を$x_t - F_t x_{t-1}$に関するものとして見なした

$$
\begin{align}
    x_t - F_t x_{t-1} | v_t &\sim \mathcal{N}(G_t v_t, \varepsilon^{2I_k}) \\
    v_t &\sim \mathcal{N}(0, Q_t) \\

    x_t - F_t x_{t-1} &\sim \mathcal{N}(0, G_t Q_t G_t^{\text{transpose}} + \varepsilon^{2I_k}) \\
    x_t|x_{t-1} &\sim \mathcal{N}(F_t x_{t-1}, G_t Q_t G_t^{\text{transpose}} + \varepsilon^{2I_k}) - \text{expr} 1 \\
        \because \text{後述の補題1より} \\
    x_{t-1}|y_{1:t-1} &\sim \mathcal{N}(x_{t-1:t-1}, V_{t-1:t-1}) \\ \text{と仮定}
\end{align}
$$

これと$\text{expr} 1$について後述の補題1を用いると

$$
\begin{align}
    \theta &= x_{t-1} \\
    \bar{\theta} &= x_{t-1:t-1} \\
    V &= V_{t-1:t-1} \\
    x &= x_t | y_{1:t-1} \\
    F &= F_t \\
    R &= G_t Q_t G_t^{\text{transpose}} + \varepsilon^{2I_k} \\

    x &\sim \mathcal{N}(F \bar{\theta}, FVF^{\text{transpose}} + R) \\
    x_t|y_{1:t-1} &\sim \mathcal{N}(F_tx_{t-1|t-1}, F_t V_{t-1:t-1} F_t^{\text{transpose}} + G_t Q_t G_t^{\text{transpose}} + \varepsilon^{2}I_k) \\

    x_t|y_{1:t-1} &\sim \mathcal{N}(x_{t|t-1}, V_{t|t-1}) - \text{expr} 2 \\
        \text{where} \\
            x_{t|t-1} &\equiv F_t x_{t-1|t-1} \\
            V_{t|t-1} &\equiv F_t V_{t-1|t-1} F_t^{\text{transpose}} + G_t Q_t G_t^{\text{transpose}}+ \epsilon^{2}I_k
\end{align}
$$

## カルマンフィルターの補題1
$$
\begin{align}
    \theta &\sim \mathcal{N}(\bar{\theta}, V) \wedge x|\theta \sim \mathcal{N}(F\theta, R) \\
    \Rightarrow x &\sim \mathcal{N}(F \bar{\theta}, FVF^{\text{transpose}} + R)
\end{align}
$$

## フィルタの導出
$$
\begin{align}
    y_t &= H_t x_t + w_t \\
        &\Rightarrow w_t = y_t - H_t x_t \\
    p(x_t|y_{1:t}) &= \\
        &\frac{p(y_t|x_t) p(x_t|y_{1:t-1})}{\int p(y_t|x_t,y_{1:t-1})p(x_t|y_{1:t-1}) dx_t}
\end{align}
$$

$p(y_t|x_t)$に着目する。\
観測モデルより
$$
\begin{align}
    y_t - H_t x_t \sim \mathcal{N}(0, R_t)
    \Rightarrow y_t|x_t \sim \mathcal{N}(H_t x_t, R_t)
\end{align}
$$
一期先予測の結果がガウス分布であると仮定する($\text{expr 2}$)
$$
    x_t|y_{1:t-1} \sim \mathcal{N}(x_{t|t-1}, V_{t|t-1})
$$

## カルマンフィルターの補題2
$$
\begin{align}
    \theta &\sim \mathcal{N}(\bar{\theta}, V) \wedge y|\theta \sim \mathcal{N}(H\theta, R) \\
    \Rightarrow \theta|y &\sim \mathcal{N}(\zeta,U) \\
    &\text{where} \\
        \zeta &= \bar{\theta} + VH^{\text{transpose}}(HVH^{\text{transpose}} + R)^{-1} (y - H \bar{\theta}) \\
        U &= V - VH^{\text{transpose}}(HVH^{\text{transpose}} + R)^{-1} HV
\end{align}
$$

補題2より
$$
\begin{align}
\theta &= x_t|y_{1:t-1} \\
\bar{\theta} &= x_{t|t-1} = F_t x_{t-1|t-1} \\
V &= V_{t|t-1} = F_t V_{t-1|t-1} F_t^{\text{transpose}} + G_t Q_t G_t^{\text{transpose}} \\
y &= y_t \\
H &= H_t \\
R &= R_t \\
    x_t|y_{1:t} &\sim \mathcal{N}(x_{t|t}, V_{t|t}) \\
    &\text{where} \\
        x_{t|t} &= x_{t|t-1} + K_t(y_t - H_t x_{t|t-1}) \\
        V_{t|t} &= V_{t|t-1} - K_t H_t V_{t|t-1} \\
        K_t &= V_{t|t-1} H_t^{\text{transpose}} (H_t V_{t|t-1} H_t^{\text{transpose}} + R_t)^{-1}
\end{align}
$$

## 補題1の証明
### 逆行列補題
$$
\begin{align}
    (V^{-1} + H^t R^{-1} H)^{-1} &= V - VH^t (R + HVH^t)^{-1} HV \quad \text{(theorem 1)} \\
    \Rightarrow (V^{-1} + H^t R^{-1} H)^{-1} H^t R^{-1} &= VH^t(R + HVH^t)^{-1} \quad \text{(theorem 2)}
\end{align}
$$

#### 証明
$$
\begin{align}
H^t R^{-1} (R + HVH^t) &= H^t + H^t R^{-1} HVH^t \quad \text{//自分の補足} \\
    &= (V^{-1} + H^t R^{-1} H)VH^t \\
H^t R^{-1} (R + HVH^t) &= (V^{-1} + H^t R^{-1} H)VH^t
\end{align}
$$
https://www.iwanttobeacat.com/entry/2019/11/24/110030

$$
\begin{align}
(I+P)^{-1} &= (I+P)^{-1} (I+P-P) \quad \text{(expr 3)} \\
    &= (I+P)^{-1} {(I+P)-P} \quad \text{(expr 4)} \\
    &= I - (I+P)^{-1} P \quad \text{(expr 5)}
\end{align}
$$

$$
\begin{align}
P + PQP &= P(I+QP) \quad \text{(expr 6)} \\
    &= (I + PQ) P \quad \text{(expr 7)}
\end{align}
$$

$$
\begin{align}
P + PQP &= P(I+QP) \\
    \Rightarrow (P+PQP)(I+QP)^{-1} &= P \\
    \Rightarrow (I+PQ)P(I+QP)^{-1} &= P \quad \text{(expr 8) true?} \\
    \Rightarrow P(I+QP)^{-1} &= (I+PQ)^{-1} P \quad \text{(expr 9) true?}
\end{align}
$$

$$
\begin{align}
P + PQP &= (I+PQ)P \\
    \Rightarrow (I+PQ)^{-1}(P+PQP) &= P \\
    \Rightarrow (I+PQ)^{-1}P(I+QP) &= P \quad \text{(expr 8) true!} \\
    \Rightarrow (I+PQ)^{-1} P &= P (I+QP)^{-1} \quad \text{(expr 9) true!}
\end{align}
$$

$$
\begin{align}
P &= (I+PQ)P(I+PQ)^{-1} \quad \text{expr 8}\longleftarrow\text{記事が間違ってる} \\
(I+PQ)^{-1}P &= P(I+PQ)^{-1} \quad \text{expr 9}\longleftarrow\text{expr 8が間違っているので当然間違っている}
\end{align}
$$

$$
\begin{align}
(A+BCD)^{-1} &= \{A(I+A^{-1}BCD)\}^{-1} \quad \text{expr 10} \\
    &= (I+A^{-1}BCD)^{-1}A^{-1} \quad \text{expr 11}\ \because\ \text{逆行列の定理 }(AB)^{-1}=B^{-1}A^{-1} \\
    &= \{I-(I+A^{-1}BCD)^{-1}A^{-1}BCD\}A^{-1} \quad \text{expr 12}\ \because\ \text{expr 5} \\
    &= A^{-1}-(I+A^{-1}BCD)^{-1}A^{-1}BCDA^{-1} \quad \text{expr 13} \\
    &= A^{-1}-A^{-1}(I+BCDA^{-1})^{-1}BCDA^{-1} \quad \text{expr 14}\ \because\ \text{expr 9 true} \\
    &= A^{-1}-A^{-1}B(I+CDA^{-1}B)^{-1}CDA^{-1} \quad \text{expr 15}\ \because\ \text{expr 9 true}
\end{align}
$$

$$
\begin{align}
\Rightarrow (I+CDA^{-1}B)^{-1}C &= (CC^{-1}+CDA^{-1}B)^{-1}C \quad \text{expr 16} \\
    &= \{C(C^{-1}+DA^{-1}B)\}^{-1}C \quad \text{expr 17} \\
    &= (C^{-1}+DA^{-1}B)^{-1}CC^{-1} \quad \text{expr 18} \\
    &= (C^{-1}+DA^{-1}B)^{-1} \quad \text{expr 19}
\end{align}
$$

$$
\therefore(A+BCD)^{-1}=A^{-1}-A^{-1}B(C^{-1}+DA^{-1}B)^{-1}DA^{-1}
$$

$$
\begin{align}
A &= V^{-1} \\
B &= H^{t} \\
C &= R^{-1} \\
D &= H \\
\text{rhs} &= V-VH^{t}(R+HVH^{t})^{-1}HV
\end{align}
$$

$$
\begin{align}
(V^{-1}+H^{t}R^{-1}H)^{-1}H^{t}R^{-1} &= VH^{t}R^{-1}-VH^{t}(R+HVH^{t})^{-1}HVH^{t}R^{-1} \\
    &= VH^{t}(R^{-1}-(R+HVH^{t})^{-1}HVH^{t}R^{-1}) \\
    &= VH^{t}(R^{-1}-R^{-1}R(R+HVH^{t}R^{-1}R)^{-1}HVH^{t}R^{-1}) \\
    &= VH^{t}(R+RR^{-1}HVH^{t})^{-1} \\
    &= VH^{t}(R+HVH^{t})^{-1}
\end{align}
$$

$$
\begin{align}
A &= R \\
B &= R \\
C &= R^{-1} \\
D &= HVH^{t} \\
\Rightarrow (V^{-1}+H^{t}R^{-1}H)^{-1} &= V-VH^{t}(R+HVH^{t})^{-1}HV
\end{align}
$$

行列式に関する定理
$$
|I_m - AB| = |I_n - BA| \quad \text{theorem 3}
$$
四分割表現したうえで単位行列とゼロ行列で埋めて恒等式とるだけ

カルマンフィルターの補題1
$$
\begin{align}
\theta &\sim \mathcal{N}(\bar{\theta}, V) \wedge x|\theta \sim \mathcal{N}(F\theta, R) \\
\Rightarrow x &\sim \mathcal{N}(F\bar{\theta}, FVF^{t} + R)
\end{align}
$$

$$
\begin{align}
p(\theta) &= \sqrt{2\pi}^{-m}|V|^{-\frac{1}{2}}\exp\left(-\tfrac{1}{2}(\theta-\bar{\theta})^{t}V^{-1}(\theta-\bar{\theta})\right) \quad \text{expr 17} \\
p(x|\theta) &= \sqrt{2\pi}^{-l}|R|^{-\frac{1}{2}}\exp\left(-\tfrac{1}{2}(x-F\theta)^{t}R^{-1}(x-F\theta)\right) \quad \text{expr 18}
\end{align}
$$

$$
p(x) = \int p(x|\theta)p(\theta)\,d\theta \quad \because\text{周辺化}
$$

$$
\begin{align}
    &p(x|\theta)p(\theta)
        = \sqrt{2\pi}^{-l-m} |V|^{-\frac{1}{2}} |R|^{-\frac{1}{2}}
        \times \text{exp}([-\frac{1}{2} (\theta-\bar{\theta})^t V^{-1} (\theta-\bar{\theta}) \\
    &\quad -\frac{1}{2} (x-F\theta)^t R^{-1} (x-F\theta)]) \quad \text{expr 19} \because \text{expr 17} \wedge \text{expr 18} \\
    &(\theta-\bar{\theta})^t V^{-1} (\theta-\bar{\theta}) + (x-F\theta)^t R^{-1} (x-F\theta) \\
    &\quad = \theta^t V^{-1} \theta - \theta^t V^{-1} \bar{\theta} - \bar{\theta}^t V^{-1} \theta + \bar{\theta}^t V^{-1} \bar{\theta} \\
    &\quad\quad + x^t R^{-1} x - x^t R^{-1} F\theta - \theta^t F^t R^{-1} x + \theta^t F^t R^{-1} F\theta \\
    &\quad = \underbrace{\theta^t(V^{-1} + F^t R^{-1} F) \theta}_{1,8項め} -\\
    &\quad\quad \underbrace{\theta^t(V^{-1} \theta + F^t R^{-1} x)}_{
        \substack{
            \theta^t(V^{-1} \bar{\theta} + F^t R^{-1} x) \\
            \text{2,7項め ここは$\bar{\theta}$が合ってる} \\
            \text{2行目2項め$\theta^t(V^{-1}\theta+\cdots$の所$\theta^t(V^{-1}\bar{\theta}+\cdots$な気がするけど俺の計算ミス？}
        }
    } -\\
    &\quad\quad \underbrace{(\bar{\theta} V^{-1} + x^t R^{-1} F)\theta}_{
        \substack{
            - (\bar{\theta}^t V^{-1} + x^t R^{-1}F)\theta\text{?} \\
            \text{3,6項め}
        }
    } +
    \underbrace{\bar{\theta}^t V^{-1} \bar{\theta} + x^t R^{-1} x}_{\text{4,5項め}} \\
    &\quad = (\theta-\zeta)^t(V^{-1} + F^t R^{-1} F)(\theta-\zeta) + \Delta \quad \text{平方完成}
\end{align}
$$

やっぱり中央2項めは$\theta^t(V^{-1}\bar{\theta}+\cdots)$だと思うし3項めは$-(\bar{\theta}^t V^{-1}\cdots)$だと思う

あと$V$と$R$は最低でも対称行列の必要がある？

where
$$
            \zeta \equiv (V^{-1} + F^t R^{-1} F)^{-1} \underbrace{(V^{-1} \bar{\theta} + F^t R^{-1} x)}_{\text{(A)}}
$$
(A): 前半の"^-1"係数ではさんでいる係数を単位行列としてキャンセルして後半の係数で項を生成している

両方$\theta$が1項目
左側$\theta$右側$\zeta$が2項目
前半部を$K^{-1}$、後半部を$L$とおくと

左側$\zeta$
$$
\begin{align}
(K^{-1} L)^t K\theta &= L^t (K^{-1})^t K\theta \\
    &= (\bar{\theta}^t (V^{-1})^t + x^t (R^{-1})^t F)(K^{-1})^t K\theta \\
    &= (\bar{\theta}^t (V^{-1})^t + x^t (R^{-1})^t F)((V^{-1})^t + F^t (R^{-1})^t F)^{-1} K\theta \\
    &= (\bar{\theta}^t (V^{-1})^t + x^t (R^{-1})^t F)((V^{-1})^t + F^t (R^{-1})^t F)^{-1} (V^{-1} + F^t R^{-1} F)\theta
\end{align}
$$
よしんば、$V$と$R$が対称行列である場合
$V = \epsilon Ik$, $R = Qt$(確率変数がそれぞれ**の成分同士が**独立の分散共分散行列なので対角行列)
具体例は共に対称行列
$$
= (\bar{\theta}^t V^{-1} + x^t R^{-1} F)\theta
$$
これを減算

両方$\zeta$
$$
\begin{align}
(K^{-1} L)^t K (K^{-1} L) &= L^t (K^{-1})^t L \\
    &= (\bar{\theta}^t (V^{-1})^t + x^t (R^{-1})^t F) (K^{-1})^t (V^{-1} \bar{\theta} + F^t R^{-1} x)
\end{align}
$$

これを加算
$$
\begin{align}
\Delta &\equiv -\zeta^t(V^{-1} + F^t R^{-1} F)\zeta + \bar{\theta}^t V^{-1} \bar{\theta} + x^t R^{-1} x \\
    &\because -L^t (K^{-1})^t K K^{-1} L = -L^t (K^{-1})^t L
\end{align}
$$
両方$\zeta$をキャンセル


残りは$\theta$も$\zeta$もかかわらないあまりの項

※ここの平方完成の手計算をすっ飛ばしています

// やりました
$$
\begin{align}
p(x) &= \sqrt{2\pi}^{-l-m} |V|^{-\frac{1}{2}} |R|^{-\frac{1}{2}} \text{exp}(-\Delta/2) \\
    &\quad \times \int \frac{\text{exp}(-(\theta-\zeta)^t(V^{-1} + F^t R^{-1} F)(\theta-\zeta)}{2} d\theta
\end{align}
$$

// ガウス積分

// $A$(ここでは$(V^{-1} + F^t R^{-1} F)$)が正定値対称な共変行列である場合
$$
\sqrt{\frac{(2\pi)^n}{det A}}
$$
// 多変数（多変量）正規分布（非退化な）
// https://manabitimes.jp/math/938

// ヤコビアン
$$
    = \sqrt{2\pi}^{-l-m} |V|^{-\frac{1}{2}} |R|^{-\frac{1}{2}} \text{exp}(-\frac{\Delta}{2}) \times \sqrt{2\pi}^m |V^{-1} + F^t R^{-1} F|^{-\frac{1}{2}}
$$

※ここの積分の手計算をすっ飛ばしています

// 手計算というより定型
$$
\begin{align}
\Delta
    &= -\left\{(V^{-1} + F^{t} R^{-1} F)^{-1}(V^{-1}\bar{\theta} + F^{t} R^{-1} x)\right\}^{t}
        (V^{-1} + F^{t} R^{-1} F)
        (V^{-1} + F^{t} R^{-1} F)^{-1}
        (V^{-1}\bar{\theta} + F^{t} R^{-1} x) \\
    &\quad + \bar{\theta}^{t} V^{-1}\bar{\theta} + x^{t} R^{-1} x \\
    &= -\left\{(V^{-1} + F^{t} R^{-1} F)^{-1}(V^{-1}\bar{\theta} + F^{t} R^{-1} x)\right\}^{t}
        (V^{-1}\bar{\theta} + F^{t} R^{-1} x) \\
    &\quad + \bar{\theta}^{t} V^{-1}\bar{\theta} + x^{t} R^{-1} x \\
    &= -(V^{-1}\bar{\theta} + F^{t} R^{-1} x)^{t}
        \left((V^{-1} + F^{t} R^{-1} F)^{-1}\right)^{t}
        (V^{-1}\bar{\theta} + F^{t} R^{-1} x) \\
    &\quad + \bar{\theta}^{t} V^{-1}\bar{\theta} + x^{t} R^{-1} x \\
    &= -(\bar{\theta}^{t} (V^{-1})^{t} + x^{t} (R^{-1})^{t} F)
        \left((V^{-1})^{t} + F^{t} (R^{-1})^{t} F\right)^{-1}
        (V^{-1}\bar{\theta} + F^{t} R^{-1} x) \\
    &\quad + \bar{\theta}^{t} V^{-1}\bar{\theta} + x^{t} R^{-1} x
\end{align}
$$
    // VとRが対称行列であると仮定
$$
\begin{align}
    &= -(\bar{\theta}^{t} V^{-1} + x^{t} R^{-1} F)
        (V^{-1} + F^{t} R^{-1} F)^{-1}
        (V^{-1}\bar{\theta} + F^{t} R^{-1} x) \\
    &\quad + \bar{\theta}^{t} V^{-1}\bar{\theta} + x^{t} R^{-1} x \\
    &= -(\bar{\theta}^{t} V^{-1})(V^{-1} + F^{t} R^{-1} F)^{-1}(V^{-1}\bar{\theta}) \\
    &\quad - (\bar{\theta}^{t} V^{-1})(V^{-1} + F^{t} R^{-1} F)^{-1}(F^{t} R^{-1})x \\
    &\quad - x^{t}(R^{-1} F)(V^{-1} + F^{t} R^{-1} F)^{-1}(V^{-1}\bar{\theta}) \\
    &\quad - x^{t}(R^{-1} F)(V^{-1} + F^{t} R^{-1} F)^{-1}(F^{t} R^{-1})x \\
    &\quad + \bar{\theta}^{t} V^{-1}\bar{\theta} + x^{t} R^{-1} x \\
    &= x^{t}\left[R^{-1} - (R^{-1} F)(V^{-1} + F^{t} R^{-1} F)^{-1}(F^{t} R^{-1})\right]x \\
    &\quad - (\bar{\theta}^{t} V^{-1})(V^{-1} + F^{t} R^{-1} F)^{-1}(F^{t} R^{-1})x \\
    &\quad - x^{t}(R^{-1} F)(V^{-1} + F^{t} R^{-1} F)^{-1}(V^{-1}\bar{\theta}) \\
    &\quad + \bar{\theta}^{t} V^{-1}\bar{\theta}
        - (\bar{\theta}^{t} V^{-1})(V^{-1} + F^{t} R^{-1} F)^{-1}(V^{-1}\bar{\theta})
\end{align}
$$
    // ここで$(V^{-1} + H^{t} R^{-1} H)^{-1} H^{t} R^{-1} = VH^{t}(R + HVH^{t})^{-1}$
    // $(V^{-1} + H^{t} R^{-1} H)^{-1} = V - VH^{t}(R + HVH^{t})^{-1} HV$を用いて
$$
\begin{align}
    &\text{theorem 2の逆} \\
    &\quad V = R^{-1},\quad H = F^{t},\quad R = V^{-1} \\
    &\text{theorem 1} \\
    &\quad V = V,\quad H = F,\quad R = R \\
    &\text{theorem 1の逆} \\
    &\quad V = R^{-1},\quad H = F^{t},\quad R = V^{-1} \\
    &\text{theorem 2} \\
    &\quad V = V,\quad H = F,\quad R = R
\end{align}
$$

$$
\begin{align}
    &= x^{t}(R + FVF^{t})^{-1} x
        - \bar{\theta}^{t} V^{-1} VF^{t}(R + FVF^{t})^{-1} x \\
    &\quad - x^{t}(R + FVF^{t})^{-1} FV V^{-1}\bar{\theta}
        - (\bar{\theta}^{t} V^{-1})\left[V - VF^{t}(R + FVF^{t})^{-1} FV\right](V^{-1}\bar{\theta}) \\
    &\quad + \bar{\theta}^{t} V^{-1}\bar{\theta} \\
    &= x^{t}(R + FVF^{t})^{-1} x
        - \bar{\theta}^{t} F^{t}(R + FVF^{t})^{-1} x \\
    &\quad - x^{t}(R + FVF^{t})^{-1} F\bar{\theta}
        - (\bar{\theta}^{t} V^{-1} V V^{-1}\bar{\theta}) \\
    &\quad + (\bar{\theta}^{t} V^{-1} VF^{t}(R + FVF^{t})^{-1} FV V^{-1}\bar{\theta})
        + \bar{\theta}^{t} V^{-1}\bar{\theta} \\
    &= x^{t}(R + FVF^{t})^{-1} x
        - \bar{\theta}^{t} F^{t}(R + FVF^{t})^{-1} x \\
    &\quad - x^{t}(R + FVF^{t})^{-1} F\bar{\theta}
        - (\bar{\theta}^{t} V^{-1}\bar{\theta}) \\
    &\quad + (\bar{\theta}^{t} F^{t}(R + FVF^{t})^{-1} F\bar{\theta})
        + \bar{\theta}^{t} V^{-1}\bar{\theta} \\
    &= x^{t}(R + FVF^{t})^{-1} x
        - \bar{\theta}^{t} F^{t}(R + FVF^{t})^{-1} x \\
    &\quad - x^{t}(R + FVF^{t})^{-1} F\bar{\theta}
        + (\bar{\theta}^{t} F^{t}(R + FVF^{t})^{-1} F\bar{\theta}) \\
    &= (x - F\bar{\theta})^{t}(FVF^{t} + R)^{-1}(x - F\bar{\theta})
\end{align}
$$
// ここまで自分の
$$
\begin{align}
    &= x^{t}(R + FVF^{t})^{-1} x
        - \bar{\theta}^{t} V^{-1} V F^{t}(FVF^{t} + R)^{-1} x \\
    &\quad - x^{t}(FVF^{t} + R)^{-1} FV V^{-1}\bar{\theta} \\
    &\quad + \bar{\theta}^{t} V^{-1}\bar{\theta}
        - \bar{\theta}^{t} V^{-1}\left[V - VF^{t}(FVF^{t} + R)^{-1} FV\right]V^{-1}\bar{\theta} \\
    &= (x - F\bar{\theta})^{t}(FVF^{t} + R)^{-1}(x - F\bar{\theta})
        \quad x\text{について平方完成}
\end{align}
$$
※ここの平方完成の手計算をすっ飛ばしています
// やりました

$$
\begin{align}
|V||R||V^{-1} + F^{t} R^{-1} F|
    &= |R||I_m + VF^{t} R^{-1} F| \\
    &= |R||I_l + R^{-1} FVF^{t}| \\
    &= |R + FVF^{t}|
\end{align}
$$
※ここの1から2行目の手計算をすっ飛ばしています
// すっ飛ばすも何もそのままです

全部代入して終わり


カルマンフィルターの補題2
$$
\begin{align}
    &\theta \sim \mathcal{N}(\bar{\theta}, V) \wedge y|\theta \sim \mathcal{N}(H\theta, R) \\
    &\Rightarrow \theta|y \sim \mathcal{N}(\zeta, U) \\
    &\text{where} \\
    &\quad \zeta = \bar{\theta} + VH^{\text{transpose}}(HVH^{\text{transpose}} + R)^{-1}(y - H\bar{\theta}) \\
    &\quad U = V - VH^{\text{transpose}}(HVH^{\text{transpose}} + R)^{-1} HV
\end{align}
$$

$$
\begin{align}
p(\theta) &= \sqrt{2\pi}^{-m} |V|^{-1/2} \exp\left(-\frac{1}{2}(\theta - \bar{\theta})^{t} V^{-1}(\theta - \bar{\theta})\right) \\
p(y|\theta) &= \sqrt{2\pi}^{-l} |R|^{-1/2} \exp\left(-\frac{1}{2}(y - H\theta)^{t} R^{-1}(y - H\theta)\right)
\end{align}
$$

$$
p(\theta|y) = \frac{p(y|\theta)p(\theta)}{\int p(y|\theta)p(\theta)\,d\theta} \quad \because \text{ベイズ}
$$

$$
\begin{align}
p(y|\theta)p(\theta)
    &= \sqrt{2\pi}^{-l-m} |V|^{-1/2} |R|^{-1/2} \\
    &\quad \times \exp\left(-\frac{1}{2}(\theta - \bar{\theta})^{t} V^{-1}(\theta - \bar{\theta})
        - \frac{1}{2}(y - H\theta)^{t} R^{-1}(y - H\theta)\right)
\end{align}
$$

$$
\begin{align}
    &(\theta - \bar{\theta})^{t} V^{-1}(\theta - \bar{\theta}) + (y - H\theta)^{t} R^{-1}(y - H\theta) \\
    &\quad = \theta^{t} V^{-1}\theta - \theta^{t} V^{-1}\bar{\theta} - \bar{\theta}^{t} V^{-1}\theta + \bar{\theta}^{t} V^{-1}\bar{\theta} \\
    &\quad\quad + y^{t} R^{-1} y - y^{t} R^{-1} H\theta - \theta^{t} H^{t} R^{-1} y + \theta^{t} H^{t} R^{-1} H\theta \\
    &\quad = \theta^{t}(V^{-1} + H^{t} R^{-1} H)\theta
        - \theta^{t}(V^{-1}\bar{\theta} + H^{t} R^{-1} y) \\
    &\quad\quad - (\bar{\theta}^{t} V^{-1} + y^{t} R^{-1} H)\theta
        + \bar{\theta}^{t} V^{-1}\bar{\theta} + y^{t} R^{-1} y \\
    &\quad = (\theta - \zeta)^{t}(V^{-1} + H^{t} R^{-1} H)(\theta - \zeta) + \Delta \\
    &\quad\quad \text{where} \\
    &\quad\quad\quad \zeta = (V^{-1} + H^{t} R^{-1} H)^{-1}(V^{-1}\bar{\theta} + H^{t} R^{-1} y) \\
    &\quad\quad\quad \Delta = -\zeta^t(V^{-1} + H^{t} R^{-1} H)\zeta + \bar{\theta}^{t} V^{-1}\bar{\theta} + y^{t} R^{-1} y
\end{align}
$$
$P(\theta|y)$に代入、ガウス積分から
$$
\begin{align}
p(\theta|y)
    &= \frac{\exp\left(-\frac{1}{2}(\theta - \zeta)^{t}(V^{-1} + H^{t} R^{-1} H)(\theta - \zeta) - \frac{1}{2}\Delta\right)}
        {\int \exp\left(-\frac{1}{2}(\theta - \zeta)^{t}(V^{-1} + H^{t} R^{-1} H)(\theta - \zeta) - \frac{1}{2}\Delta\right) d\theta} \\
    &= \frac{\exp\left(-\frac{1}{2}(\theta - \zeta)^{t}(V^{-1} + H^{t} R^{-1} H)(\theta - \zeta)\right)}
        {\sqrt{2\pi}^{m} \left|(V^{-1} + H^{t} R^{-1} H)\right|^{-\frac{1}{2}}} \\
    &= \sqrt{2\pi}^{-m} |U|^{-\frac{1}{2}} \exp\left(-\frac{1}{2}(\theta - \zeta)^{t} U^{-1}(\theta - \zeta)\right) \\
    &\quad\quad \text{where} \\
    &\quad\quad\quad U \equiv (V^{-1} + H^{t} R^{-1} H)^{-1}
\end{align}
$$

$$
\begin{align}
\zeta
    &= (V^{-1} + H^{t} R^{-1} H)^{-1}(V^{-1}\bar{\theta} + H^{t} R^{-1} y) \\
    &= (V^{-1} + H^{t} R^{-1} H)^{-1} V^{-1}\bar{\theta}
        + (V^{-1} + H^{t} R^{-1} H)^{-1} H^{t} R^{-1} y
\end{align}
$$

$$
\begin{align}
    &\text{theorem 2} \\
    &\quad V = V,\quad H = H,\quad R = R \\
    &\text{theorem 1} \\
    &\quad V = V,\quad H = H,\quad R = R
\end{align}
$$

$$
\begin{align}
    &= (V - VH^{t}(R + HVH^{t})^{-1} HV)V^{-1}\bar{\theta}
        + VH^{t}(R + HVH^{t})^{-1} y \\
    &= \bar{\theta} - VH^{t}(R + HVH^{t})^{-1} H\bar{\theta}
        + VH^{t}(R + HVH^{t})^{-1} y \\
    &= \bar{\theta} + VH^{t}(R + HVH^{t})^{-1}(y - H\bar{\theta})
\end{align}
$$

$$
\begin{align}
U
    &= (V^{-1} + H^{t} R^{-1} H)^{-1} \\
    &\quad\quad \text{theorem 2} \\
    &\quad\quad\quad V = V,\quad H = H,\quad R = R \\
    &= V - VH^{t}(R + HVH^{t})^{-1} HV
\end{align}
$$

終わり

ようやっと付録の証明も含めて全部手で追えた、おつです