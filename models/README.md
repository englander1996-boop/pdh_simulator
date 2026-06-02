# models/ — 学習済みサロゲートモデル

## 役割

蒸留塔 Dist1 / Dist3 を HYSYS の代わりに高速予測するための **学習済み GPR (ガウス過程回帰)
サロゲートモデル**。HYSYS 解とほぼ一致するよう学習されており、HYSYS COM を起動せずに
HYSYS 精度・FUG 並みの速度で塔を評価できる。`main.py` / `sub2.py` / `exp3.py` 等の
SM バックエンドに必須。

## ファイル一覧

| ファイル | 内容 |
|---|---|
| `column1_sm.pkl` | Dist1 (塔1) の GPR サロゲート (dict)。学習域 N_stages = 30〜60。 |
| `column3_sm.pkl` | Dist3 (塔3, C3 スプリッタ) の GPR サロゲート (dict)。 |

## 使い方・位置づけ

- `src/distillation_sm.py` がロードする (`solve_column1_via_sm` / `solve_column3_via_sm`)。
  プロセス内キャッシュ (`_MODEL_CACHE`) で 1 回だけロードされる。
- これらが無いと SM バックエンド (Dist1=SM / Dist3=SM) での評価ができない。Dist2 は SM 化できず
  HYSYS COM (`main.py`) または in-house rigorous (`sub2.py`) で評価する。
- 詳細は `src/SPEC_distillation.md` を参照。
