import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import glob
import curv_fit as cf
import sys

def get_step(target, target_library):
    # 1から9までのstabilityのlistを返す
    try:
        target = target + "/"
        return cf.stability_divided_by_base(target, target_library)
    except Exception as e:
        print(f"Error in get_step: {e}")
        return []

def plot(data, target_library, temp, color):
    # DataFrameを作成
    stabilities_df = pd.DataFrame.from_dict(data[temp], orient='index', columns=['stability'])
    stabilities_df.index = pd.MultiIndex.from_tuples(stabilities_df.index, names=['step', 'target'])
    
    # インデックスをリセットして、カラムに分ける
    stabilities_df.reset_index(inplace=True)
    
    # 'temp' カラムを追加
    stabilities_df['temp'] = temp

    print(stabilities_df.head())  # 確認のため最初の数行を表示
    print(stabilities_df.dtypes)  # 確認のためのデータ型の表示

    plt.title(f'Stability Analysis for {target_library}')
    
    # hue='temp' を維持し、カラーパレットを指定
    sns.lineplot(data=stabilities_df, x='step', y='stability', hue='temp', palette=[color])
    
    plt.show()

def main():
    if len(sys.argv) != 4:
        print("usage: python3 comparison_curv_fit.py <type of l> <target>")
    steps = list(range(1, 10))
    target_library = sys.argv[1]
    target = sys.argv[2]
    targets = glob.glob(f"/home/user/SA-EDS/{target}/{target_library}_*/*")
    # targets = glob.glob(f"/home/user/SA-EDS/{target}/{target_library}_initial_*/{target_library}-*")
    print(f"/home/user/SA-EDS/{target}/{target_library}_*/*")
    # print(f"/home/user/SA-EDS/{target}/{target_library}_initial_*/{target_library}-*")
    print(f"Targets found: {targets}")
    
    # データを作成
    data = {}
    
    for target in targets:
        temp = target.split("_")[-2]
        print(temp)

        stabilities = get_step(target, target_library)
        if len(stabilities) != 9:
            continue
        print(f"Stabilities for {target}: {stabilities}")
        for step in steps:
            if len(stabilities) >= step:
                if temp not in data:
                    data[temp] = {}
                data[temp][(step, target.split("/")[-1])] = stabilities[step-1]  

    print(f"Data dictionary: {data}")

    if not data:
        print("No data to plot.")
        return

    # カスタムカラーリストを作成

    temp_list = ["277", "298", "308", "318", "328", "338"]
    temp_list = ["277", "308", "328"]
    if sys.argv[2] == "int_second":
       temp_list = ["277", "308", "328"]
    color_list = [(0.9692676429840149, 0.3347286967411282, 0.6175245379502425), \
                (0.8826628895184139, 0.6205292536652166, 0.13305832188134544), \
                (0.45778524885227995, 0.8205005154203934, 0.28102349811240404), \
                (0.14330890377231592, 0.7961171141278752, 0.6749754171147194), \
                (0.15561284657119115, 0.6357382064596534, 0.927372835549898), \
                (0.5757357436608456, 0.43914632757667336, 0.8581768581768582)]
  # 'husl' で6色を生成
    for i, temp in enumerate(temp_list):
        plot(data, target_library, temp, color_list[i])  # 色をリストから取得
    
    plt.savefig(f'/home/user/SA-EDS/fig/curv_{target_library}_{sys.argv[2]}_3.png') 
    print(f'/home/user/SA-EDS/fig/curv_{target_library}_{sys.argv[2]}_3.png')

if __name__ == '__main__':
    main()
