## Load python packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import GridSearchCV
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc, plot_roc_curve


## Load and pre-process data
data = pd.read_csv('BreastCancer.csv', sep=',')
feature = [
    'Cl.thickness', 'Cell.size', 'Cell.shape', 'Marg.adhesion', 'Epith.c.size',
    'Bare.nuclei', 'Bl.cromatin', 'Normal.nucleoli', 'Mitoses'
]
X = data[feature]
Y = np.array(data['Class'].replace('benign', 0).replace('malignant', 1))

# 利用均值进行空缺值填充
X_mean = X.mean(axis=0)
for i in range(len(X.T)):
    X.iloc[:, i] = X.iloc[:, i].fillna(X_mean[i])

# 使用standard/z-score scaling 对数据做scaling
X_sd = StandardScaler().fit_transform(X)
X = pd.DataFrame(X_sd)
X.columns = feature

## Dividing data
# Whole Data Set → Discovery Set/Validation Set
random_state = np.random.RandomState(1289237)
X_discovery, X_validation, y_discovery, y_validation = train_test_split(
    X, Y, test_size=0.2, random_state=random_state)
X_discovery.index = range(len(X_discovery))
X_validation.index = range(len(X_validation))
print('number of discovery samples: {}, validation samples: {}'.format(
    X_discovery.shape[0], X_validation.shape[0]))

# Discovery set → Training Set/Test Set
# Using K-fold Cross-validation


## Feature selection based on cross-validation
# 模型函数(clf_select)
def clf_select(name):
    if name == 'DT':
        clf = DecisionTreeClassifier(max_depth=5,
                                     min_samples_leaf=5,
                                     criterion='gini')
    elif name == 'DT_cv':
        tree_para = {'max_depth': [3, 5, 7, 9]}
        clf = GridSearchCV(DecisionTreeClassifier(), tree_para, cv=5)
    elif name == 'SVM':
        clf = SVC(kernel='rbf', probability=True, C=1)
    elif name == 'SVM_cv':
        tree_para = {
            'C': [
                0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000,
                100000
            ]
        }
        clf = GridSearchCV(SVC(kernel='rbf', probability=True),
                           tree_para,
                           cv=5)
    elif name == 'RF':
        clf = RandomForestClassifier(n_estimators=50, max_depth=5)
    elif name == 'RF_cv':
        tree_para = {'n_estimators': [25, 50, 75], 'max_depth': [3, 4, 5]}
        clf = GridSearchCV(RandomForestClassifier(), tree_para, cv=5)
    elif name == 'LR_cv':
        tree_para = {
            'C': [
                0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000,
                100000
            ]
        }
        clf = GridSearchCV(LogisticRegression(penalty='l2',
                                              solver='liblinear'),
                           tree_para,
                           cv=5)
    return clf


# 计算可能的特征组合
from itertools import combinations
feature_list = []
for i in combinations(feature, 3):
    feature_list.append(list(i))

# 对每个特征组合进行交叉验证计算5个子Test Sets的平均AUC
result_train = pd.DataFrame(columns={'feature', 'AUC_mean'})
result_test = pd.DataFrame(columns={'feature', 'AUC_mean'})

for j in range(len(feature_list)):
    #print(j)
    # 交叉验证划分
    skf = StratifiedKFold(n_splits=5, random_state=1, shuffle=True)
    # 交叉验证中每一折结果
    result_train_ = pd.DataFrame(columns={'num', 'AUC'})
    result_test_ = pd.DataFrame(columns={'num', 'AUC'})
    n = 0
    for train, test in skf.split(list(X_discovery.index), y_discovery):
        # Training and Test Set
        X_train = X_discovery.loc[train, feature_list[j]]
        X_test = X_discovery.loc[test, feature_list[j]]
        y_train = y_discovery[train]
        y_test = y_discovery[test]

        # 模型训练，我们自定义clf_select函数
        clf = clf_select('DT_cv')
        clf.fit(X_train, y_train)

        # 模型预测结果
        pred_proba_train = clf.predict_proba(X_train)
        fpr_train, tpr_train, thresholds = roc_curve(y_train,
                                                     pred_proba_train[:, 1])
        roc_auc_train = auc(fpr_train, tpr_train)
        pred_proba_test = clf.predict_proba(X_test)
        fpr_test, tpr_test, thresholds = roc_curve(y_test, pred_proba_test[:,
                                                                           1])
        roc_auc_test = auc(fpr_test, tpr_test)

        result_train_.loc[n, 'num'] = n
        result_train_.loc[n, 'AUC'] = roc_auc_train
        result_test_.loc[n, 'num'] = n
        result_test_.loc[n, 'AUC'] = roc_auc_test
        n = n + 1

    # 模型Test Set平均AUC计算
    result_train.loc[j, 'feature'] = ','.join(feature_list[j])
    result_train.loc[j, 'AUC_mean'] = result_train_['AUC'].mean()
    result_test.loc[j, 'feature'] = ','.join(feature_list[j])
    result_test.loc[j, 'AUC_mean'] = result_test_['AUC'].mean()

# 根据5个子Test Sets的平均AUC选择最佳特征组合
best_feature = result_test.loc[
    result_test.sort_values('AUC_mean', ascending=False).index[0],
    'feature'].split(',')
print(best_feature)

## Evaluate on validation set
# 我们首先利用cross-validation选择好的超参数重新对整个Discovery set进行训练得到一个预测模型
clf = clf_select('DT_cv')
clf.fit(X_discovery.loc[:, best_feature], y_discovery)

# 绘制ROC
plot_roc_curve(clf, X_validation.loc[:, best_feature], y_validation)
plt.show()
