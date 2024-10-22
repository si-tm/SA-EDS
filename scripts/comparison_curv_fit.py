import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import glob
import curv_fit as cf

def get_step(target, target_library):
    # 1から9までのstabilityのlistを返す
    try:
        target = target + "/"
        return cf.stability(target, target_library)
    except Exception as e:
        print(f"Error in get_step: {e}")
        return []

def main():
    steps = list(range(1, 10))
    target_library = "L1"
    targets = glob.glob(f"/home/user/SA-EDS/int_initial/{target_library}_initial_*/{target_library}-*")
    # targets = glob.glob(f"/home/user/SA-EDS/int_initial/{target_library}_initial_*/{target_library}-GA100000-0.80-MSS-9_277*")

    print(f"Targets found: {targets}")
    
    # 乗客数のデータを作成（適当なサンプルデータを使用）
    data = {}
    
    for target in targets:
        stabilities = get_step(target, target_library)
        if len(stabilities) != 9:
            continue
        print(f"Stabilities for {target}: {stabilities}")
        for step in steps:
            if len(stabilities) >= step:
                data[(step, target.split("/")[-1])] = stabilities[step-1]  

    print(f"Data dictionary: {data}")
    


    if not data:
        print("No data to plot.")
        return

    # DataFrameを作成
    stabilities_df = pd.DataFrame.from_dict(data, orient='index', columns=['stability'])
    stabilities_df.index = pd.MultiIndex.from_tuples(stabilities_df.index, names=['step', 'target'])
    
    # インデックスをリセットして、カラムに分ける
    stabilities_df.reset_index(inplace=True)

    print(stabilities_df.head())  # 確認のため最初の数行を表示
    print(stabilities_df.dtypes)  # 確認のため最初の数行を表示

    plt.title(f'Stability Analysis for {target_library}')
    
    # プロット
    sns.lineplot(data=stabilities_df, x='step', y='stability')
    # sns.lineplot(data=stabilities_df, x='step', y='stability', hue='target')
    plt.show()
    plt.savefig(f'/home/user/SA-EDS/fig/curv_{target_library}.png') 




if __name__ == '__main__':
    main()
