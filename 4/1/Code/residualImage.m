function output = residualImage(corruptedData, biasRemoved, biasField)
output = corruptedData - (biasRemoved .* biasField);