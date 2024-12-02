import os
import pickle  # モデルロードに使用
import graphviz  # Graphvizをインポート
from sklearn import tree  # 決定木の可視化に必要

def load_model(path):
    """モデルを指定パスからロード"""
    with open(path, 'rb') as f:
        model = pickle.load(f)
    return model

if __name__ == '__main__':
    # モデルをロード
    model_path = "../saved_model/bagging_model_L1_int_initial.pkl"
    bagging_model = load_model(model_path)

    # Bagging モデルから最初の推定器 (RandomForest) を取得
    base_estimators = bagging_model.estimators_  # Baggingモデル内のベース推定器
    random_forest = base_estimators[0]  # 最初の RandomForest を取得

    # RandomForest 内の最初の決定木を取得
    first_tree = random_forest.estimators_[0]  # RandomForest内の最初の決定木

    # Graphviz形式でエクスポート
    dot_data = tree.export_graphviz(
        first_tree,
        out_file=None,  # ファイルではなく文字列として出力
        filled=True,  # ノードを色付け
        rounded=True,  # ノードを丸みを帯びた形に
        special_characters=True  # 特殊文字の扱い
    )

    # グラフを画像ファイルとして保存
    graph = graphviz.Source(dot_data)
    output_dir = '../fig'  # 保存ディレクトリ
    os.makedirs(output_dir, exist_ok=True)  # ディレクトリを作成
    graph.format = 'png'
    graph.render(os.path.join(output_dir, 'decision_tree'))  # 保存
