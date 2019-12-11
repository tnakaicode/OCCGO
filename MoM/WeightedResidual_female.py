import matplotlib.pyplot as plt
import pandas
from pandas import Series
from pandas import DataFrame
from pandas import concat
from statsmodels.tsa.ar_model import AR
from sklearn.metrics import mean_squared_error

series = pandas.Series.from_csv('daily-total-female-births.csv', header=0)
print(series.head())

# create lagged dataset
values = DataFrame(series.values)
dataframe = concat([values.shift(1), values], axis=1)
dataframe.columns = ['t-1', 't+1']

# split into train and test sets
X = dataframe.values
train_size = int(len(X) * 0.66)
train, test = X[1:train_size], X[train_size:]
train_X, train_y = train[:, 0], train[:, 1]
test_X, test_y = test[:, 0], test[:, 1]

# persistence model
predictions = [x for x in test_X]

# calculate residuals
residuals = [test_y[i] - predictions[i] for i in range(len(predictions))]
residuals = DataFrame(residuals)
print(residuals.head())

# persistence model on training set
train_pred = [x for x in train_X]

# calculate residuals
train_resid = [train_y[i] - train_pred[i] for i in range(len(train_pred))]

# model the training set residuals
model = AR(train_resid)
model_fit = model.fit()
window = model_fit.k_ar
coef = model_fit.params
print('Lag=%d, Coef=%s' % (window, coef))

# walk forward over time steps in test
history = train_resid[len(train_resid) - window:]
history = [history[i] for i in range(len(history))]
predictions = list()
expected_error = list()
for t in range(len(test_y)):
    # persistence
    yhat = test_X[t]
    error = test_y[t] - yhat
    expected_error.append(error)
    
    # predict error
    length = len(history)
    lag = [history[i] for i in range(length - window, length)]
    pred_error = coef[0]
    for d in range(window):
        pred_error += coef[d + 1] * lag[window - d - 1]

    # correct the prediction
    yhat = yhat + pred_error
    predictions.append(yhat)
    history.append(error)
    print('predicted=%f, expected=%f' % (yhat, test_y[t]))

# error
mse = mean_squared_error(test_y, predictions)
print('Test MSE: %.3f' % mse)

# plot predicted error
plt.figure()
plt.plot(test_y)
plt.plot(predictions, color='red')
plt.show()
