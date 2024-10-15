#!/bin/bash

# 使用率の閾値を指定 (例: 20%以下の時に実行する)
THRESHOLD=20

# 実行したいコマンド (例: echo コマンド)
COMMAND="echo 'CPU usage is below threshold. Running command...'"

# CPU使用率を取得する関数
get_cpu_usage() {
    # 1回目のCPU統計を取得
    CPU=($(cat /proc/stat | grep '^cpu ')) # 全CPUの統計を取得
    IDLE_OLD=${CPU[4]} # アイドル時間
    TOTAL_OLD=0
    for VALUE in "${CPU[@]:1}"; do
        TOTAL_OLD=$((TOTAL_OLD + VALUE))
    done

    # 1秒待ってから2回目のCPU統計を取得
    sleep 1
    CPU=($(cat /proc/stat | grep '^cpu '))
    IDLE_NEW=${CPU[4]}
    TOTAL_NEW=0
    for VALUE in "${CPU[@]:1}"; do
        TOTAL_NEW=$((TOTAL_NEW + VALUE))
    done

    # 使用率を計算
    DIFF_IDLE=$((IDLE_NEW - IDLE_OLD))
    DIFF_TOTAL=$((TOTAL_NEW - TOTAL_OLD))
    CPU_USAGE=$((100 * (DIFF_TOTAL - DIFF_IDLE) / DIFF_TOTAL))

    echo $CPU_USAGE
}

while true; do
    # CPU使用率を取得
    CPU_USAGE=$(get_cpu_usage)

    # CPU使用率が閾値より低いかチェック
    if (( CPU_USAGE < THRESHOLD )); then
        # CPU使用率が閾値より低ければコマンドを実行
        eval $COMMAND
        break
    fi

    # 1秒待つ
    sleep 1
done
