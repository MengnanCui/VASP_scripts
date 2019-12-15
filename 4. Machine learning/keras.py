#!/gpfs/home/mncui/soft/anaconda3/bin/python
import numpy as np
import pandas as pd
import keras
from keras.models import Sequential
from keras.layers import Dense
from keras.layers import CuDNNLSTM
from keras.wrappers.scikit_learn import KerasRegressor
from sklearn.cross_validation import cross_val_score
from sklearn.cross_validation import KFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
import matplotlib.pyplot as plt
from keras.layers import Dropout
log = open('run.log','a')


calina = pd.read_table('symple.dat',header = None)
x = calina.values[0:4000000,0:50]
#calina = calina.values[:,0:25]
#calinb = pd.rad_table('../pos10/symple.dat',header = None)
#calinb = calinb.values[:,:]
#calin = np.r_[calina,calinb]
caloua = pd.read_csv('chgcar.dat',header = None)
y = calina.values[0:4000000,0]
x_val = calina.values[4000001:,0:50]
y_val = caloua.values[4000001:,0]
#caloub = pd.rad_csv('../pos10/chgcar.dat',header = None)
#caloub = calinb.values[:,:]
#calou = np.r_[caloua,caloub]

model = Sequential()
model.add(Dropout(0.5, input_shape=(50,)))
model.add(Dense(512,input_dim=50, init='normal', activation='relu'))
model.add(Dropout(0.5))
model.add(Dense(512, init='normal', activation='relu'))
model.add(Dropout(0.5))
model.add(Dense(512, init='normal', activation='relu'))
model.add(Dropout(0.5))
model.add(Dense(512, init='normal', activation='relu'))
model.add(Dropout(0.5))

model.compile(loss='mean_squared_error',
              optimizer='adam',metrics=['accuracy'])
#model.summary()

#                   validation_split=0.9,
history = model.fit(x,y,epochs=100,batch_size=200000,verbose=1,validation_data=(x_val,y_val))
print(history.history.keys())
loss = history.history['loss']
val_loss = history.history['val_loss']


epochs = range(1, len(loss)+1)

plt.plot(epochs, loss, 'bo', label='Training loss')
plt.plot(epochs, val_loss, 'b', label='Validation loss')
plt.title('Training and validation loss')
plt.xlabel('Epochs')
plt.ylabel('Loss')
plt.legend()
plt.savefig('loss.png')
plt.clf()

accuracy = history.history['acc']
val_acc = history.history['val_acc']

plt.plot(epochs, accuracy,'bo', label='Accuracy')
plt.plot(epochs, val_acc,'b', label='Validation accuracy')
plt.title('Training and validation accuracy')
plt.xlabel('epochs')
plt.ylabel('accuracy')
plt.legend()
plt.savefig('accuracy.png')
