require 'torch'
require 'hdf5'


print 'Loading dataset'

train_file1 = '7_train.fasta.ref.h5'
valid_file1 = '7_valid.fasta.ref.h5'
train_file2 = '7_train.label.ref.h5'
valid_file2 = '7_valid.label.ref.h5'


loaded1 = hdf5.open(train_file1,'r')
loaded2 = hdf5.open(train_file2,'r')
tr_size = (#loaded2:read('/testdata'):all():transpose(2,1))[1]
trainData = {
    data = loaded1:read('/testxdata'):all(),
    labels = loaded2:read('/testdata'):all():transpose(2,1),
    size = function() return tr_size end
}

loaded1 = hdf5.open(valid_file1,'r')
loaded2 = hdf5.open(valid_file2,'r')
te_size = (#loaded2:read('/testdata'):all():transpose(2,1))[1]
validData = {
    data = loaded1:read('/testxdata'):all(),
    labels = loaded2:read('/testdata'):all():transpose(2,1),
    size = function() return te_size end
}
   
noutputs = (#validData.labels)[2]

print 'Finished loading dataset'






