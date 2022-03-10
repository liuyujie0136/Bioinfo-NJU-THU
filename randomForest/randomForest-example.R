library(randomForest)
library(ROCR)

data(iris)

# 去除类别为setosa的样本
df <- iris[iris$Species != 'setosa', ]
# 去除最后一列Species，产生输入矩阵
all_data <- df[, 1:(ncol(df) - 1)]
# 需要预测的类别向量。factor()函数用于把原来的3种类别编程2中类别
all_classes <- factor(df$Species)

# Calculate the size of training sets
n_samples <- nrow(all_data)
n_train <- floor(n_samples * 0.8)  

# 产生所有样本序号的重新排列
set.seed(1234)
indices <- sample(1:n_samples)
# 从以上随机排列中选取n_train个作为训练集样本序号
indices <- indices[1:n_train]

# 选取Training Set
# 以下代码从矩阵all_data按照indices里的顺序抽取相应的行拼接成一个新的矩阵
train_data <- all_data[indices,]
# 以下代码从向量all_classes中按照indices里的顺序抽取相应元素组成一个新的向量
train_classes <- all_classes[indices]
# 选取Validation Set
# -indices表示选取序号不在indices中的样本
validation_data <- all_data[-indices,]
validation_classes <- all_classes[-indices]

# training
rf_classifier = randomForest(train_data, train_classes, trees = 100)

# validation
predicted_classes <- predict(rf_classifier, validation_data)
predicted_classes
validation_classes

# confusion matrix
# 定义versicolor为正类别
positive_class <- 'versicolor'
# true positive count (TP)
TP <- sum((predicted_classes == positive_class) & (validation_classes == positive_class))
# false positive count (FP)
FP <- sum((predicted_classes == positive_class) & (validation_classes != positive_class))
# false negative count (FN)
FN <- sum((predicted_classes != positive_class) & (validation_classes == positive_class))
# true negative count (TN)
TN <- sum((predicted_classes != positive_class) & (validation_classes != positive_class))
# 构建2x2矩阵，填充以上计算的四个数
confusion <- matrix(c(TP, FP, FN, TN), nrow=2, ncol=2, byrow=TRUE)
colnames(confusion) <- c('True versicolor', 'True virginica')
rownames(confusion) <- c('Predicted versicolor', 'Predicted virginica')
confusion

# 
cat('accuracy =', (TP + TN)/(TP + TN + FP + FN), '\n')
cat('sensitivity =', TP/(TP + FN), '\n')
cat('positive predicted value =', TP/(TP + FP), '\n')
cat('specificity =', TN/(TN + FP), '\n')

# Plot ROC on validation set
predicted_probs <- predict(rf_classifier, validation_data, type = 'prob')

# 创建一个长度与测试集大小相同的0-1向量，1代表要预测的类别
validation_labels <- vector('integer', length(validation_classes))
validation_labels[validation_classes != positive_class] <- 0
validation_labels[validation_classes == positive_class] <- 1
# 通过prediction函数，使用预测为正样本的概率和真实类别创建一个对象pred
pred <- prediction(predicted_probs[, positive_class], validation_labels)

roc <- performance(pred, 'tpr', 'fpr') 
plot(roc, main = 'ROC Curve')

auc <- performance(pred, 'auc')
cat('auc =', auc@y.values[[1]], '\n')

# 特征之间相关性
GGally::ggpairs(df, columns = 1:4, ggplot2::aes(color = Species))

