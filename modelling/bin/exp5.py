import re
import os

# ==========================================
# 設定
# ==========================================
# ステップ1で保存した未計算のXMLファイルを指定
INPUT_XML = r"Z:\pdh_simulator\project\wakeup_hysys\C3C4_Splitter_Reset.xml"
# 出力するテスト用XMLファイル
OUTPUT_XML = r"Z:\pdh_simulator\project\wakeup_hysys\Test_S60_F30.xml"

# 現在のXMLに記録されている設定値（元のHYSYSでの数値）
OLD_STAGES = 40
OLD_FEED = 20

# 今回テストで生成したい新しい設定値
NEW_STAGES = 60
NEW_FEED = 30

def inject_hysys_xml():
    if not os.path.exists(INPUT_XML):
        print(f"[エラー] 入力ファイルが見つかりません: {INPUT_XML}")
        return

    print("1. XMLファイルを読み込んでいます...")
    with open(INPUT_XML, "r", encoding="utf-8", errors="ignore") as f:
        content = f.read()

    print("2. 構造的変数（トポロジー）の書き換えを実行中...")
    
    # 段数 (NumTrays または NumberOfStages) の書き換え
    # 例: <NumTrays>40</NumTrays> -> <NumTrays>60</NumTrays>
    content, stage_subs = re.subn(
        rf"(<(?:NumTrays|NumberOfStages)>){OLD_STAGES}(</(?:NumTrays|NumberOfStages)>)",
        rf"\g<1>{NEW_STAGES}\g<2>",
        content
    )
    print(f" -> 段数タグを {stage_subs} 箇所書き換えました。")

    # フィード段 (StepNumber, Stage, Location等) の書き換え
    # 例: <StepNumber>20</StepNumber> -> <StepNumber>30</StepNumber>
    content, feed_subs = re.subn(
        rf"(<(?:StepNumber|Stage|Location)>){OLD_FEED}(</(?:StepNumber|Stage|Location)>)",
        rf"\g<1>{NEW_FEED}\g<2>",
        content
    )
    print(f" -> フィード段タグを {feed_subs} 箇所書き換えました。")

    print("3. 新しいXMLファイルを保存しています...")
    with open(OUTPUT_XML, "w", encoding="utf-8") as f:
        f.write(content)
        
    print(f"\n★ XMLインジェクション完了: '{os.path.basename(OUTPUT_XML)}' を生成しました。")
    print("HYSYSでこのXMLファイルを読み込めるかテストしてください。")

if __name__ == "__main__":
    inject_hysys_xml()