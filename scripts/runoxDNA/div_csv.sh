#!/bin/bash

# 元のファイル
input_file="dir_L3_1.csv"

# 出力ファイルを初期化 (dir_L3_{0..7}.csv の内容を空にする)
for i in {0..7}; do
    > "dir_L3_1_${i}.csv"
done

# 元のファイルを1行ずつ読み込む
line_number=0
while read -r line; do
    # 行番号に対して x % 8 を計算
    i=$((line_number % 8))

    # 計算結果に基づいて対応するファイルに行を追加
    echo "$line" >> "dir_L3_1_${i}.csv"

    # 行番号をインクリメント
    line_number=$((line_number + 1))
done < "$input_file"

echo "ファイルを dir_L3_1_{0..7}.csv に分割しました。"


