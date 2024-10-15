#!/bin/bash

# ホスト名の最後の文字を取得
hostname=$(uname -n)
last_char="${hostname: -1}"

# 閾値の設定
LOAD_THRESHOLD=2.0  # 負荷平均のしきい値
CPU_IDLE_THRESHOLD=20.0  # CPUアイドル時間のしきい値（%）
MEM_FREE_THRESHOLD=2000000  # メモリの空き容量しきい値（KiB）

check_system_load() {
    # Load Averageの取得
    load_avg=$(top -bn1 | grep "load average" | awk '{print $10}' | sed 's/,//')
    
    # CPU使用率の取得
    cpu_idle=$(top -bn1 | grep "%Cpu(s)" | awk '{print $8}')
    
    # メモリ使用量の取得
    mem_free=$(top -bn1 | grep "KiB Mem" | awk '{print $4}')

    # 現在の負荷状況を表示
    echo "時刻: $(date)"
    echo "現在の負荷平均: $load_avg / しきい値: $LOAD_THRESHOLD"
    echo "現在のCPUアイドル時間: $cpu_idle% / しきい値: $CPU_IDLE_THRESHOLD%"
    echo "現在のメモリ空き容量: $mem_free KiB / しきい値: $MEM_FREE_THRESHOLD KiB"
    
    # 負荷平均がしきい値を超えているか判定
    if (( $(echo "$load_avg > $LOAD_THRESHOLD" | bc -l) )); then
        echo "高負荷: 負荷平均がしきい値を超えています。現在の負荷平均: $load_avg"
        return 1
    fi

    # CPUアイドル時間がしきい値を下回っているか判定
    if (( $(echo "$cpu_idle < $CPU_IDLE_THRESHOLD" | bc -l) )); then
        echo "高負荷: CPUアイドル時間がしきい値を下回っています。現在のアイドル時間: $cpu_idle%"
        return 1
    fi

    # メモリの空き容量がしきい値を下回っているか判定
    if [ "$mem_free" -lt "$MEM_FREE_THRESHOLD" ]; then
        echo "高負荷: メモリの空き容量がしきい値を下回っています。現在の空き容量: $mem_free KiB"
        return 1
    fi

    return 0  # 負荷は許容範囲内
}


# ファイルを読み込み、処理を開始
while read row; do
  dir=$(echo "${row}" | cut -d , -f 1)  # カンマ区切りで1列目を取得

  for i in {1..10}; do
    # フラグ変数を初期化
    command_completed=0

    while [[ $command_completed -eq 0 ]]; do
        # システム負荷をチェック
        check_system_load
        if [[ $? -eq 0 ]]; then  # 関数の戻り値を確認
            command_completed=1
            # カレントディレクトリを変更し、nohupでDockerをバックグラウンド実行
            cd /clusterhome/mhyakuzuka/oxdna/L3_1/${dir}/${dir}_${i} && nohup docker run --rm -v "$(pwd)/test:/code/oxDNA/results" oxdna > output.log 2>&1 &

            # 実行したプロセスの PID を保存
            pid=$!
            echo "Docker コンテナが PID $pid で実行中です。"

            sleep 3600  # 1時間待機
            break
        else
            echo "High load detected. Waiting for 10 minutes before retrying..."
            sleep 600  # 10分待機
        fi
    done
  done
done < "dir_L3_1_${last_char}.csv"

# すべてのバックグラウンドプロセスが終了するまで待機
wait
echo "すべての Docker コンテナの実行が完了しました。"
