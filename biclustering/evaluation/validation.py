#folder handing
import os

#data handling
import numpy as np
from scipy.stats import rankdata


## Biclustering Quality Measures

def MSR(bic):
	"""
	source:
	Cheng, Y., & Church, G. M. (2000, August). Biclustering of expression data. In Ismb (Vol. 8, No. 2000, pp. 93-103).

	Parameters:
		:param bic: numpy 2D array

	Output:
		:return: Mean Squared Residue Score
	"""
	if bic.shape[0]<=1 or bic.shape[1]<=1:
		return 0.0

	col_mean = np.mean(bic, axis=0)
	row_mean = np.mean(bic, axis=1)
	bic_mean = np.mean(bic.flatten())
	
	msr = 0

	for itr_x in range(bic.shape[0]):
		for itr_y in range(bic.shape[1]):
			msr+= (bic[itr_x,itr_y]-row_mean[itr_x]-col_mean[itr_y]+bic_mean)**2
		
	return msr/(bic.shape[0]*bic.shape[1])

def SMSR(bic):
	"""
	source:
	Mukhopadhyay, A., Maulik, U., & Bandyopadhyay, S. (2009). A novel coherence measure for discovering scaling biclusters from gene expression data. Journal of Bioinformatics and Computational Biology, 7(05), 853-868.

	Parameters:
		:param bic: numpy 2D array

	Output:
		:return: Scaled Mean Squared Residue Score
	"""
	if bic.shape[0]<=1 or bic.shape[1]<=1:
		return 0.0

	smsr = 0

	for itr_x in range(bic.shape[0]):
		for itr_y in range(bic.shape[1]):
			smsr+= ((row_mean[itr_x]*col_mean[itr_j]- bic[itr_x,itr_y]*bic_mean)**2)/(row_mean[itr_x]**2 * col_mean[itr_y]**2)
		
	return smsr/(bic.shape[0]*bic.shape[1])

def VE(bic):
	"""
	source:
	Divina, F., Pontes, B., Giráldez, R., & Aguilar-Ruiz, J. S. (2012). An effective measure for assessing the quality of biclusters. Computers in biology and medicine, 42(2), 245-256
    Parameters:
		:param bic: numpy 2D array

	Output:
		:return : Virtual Error Score
	"""
	if bic.shape[0]<=1 or bic.shape[1]<=1:
		return 0.0

	rho = np.mean(bic,axis=0)
	rho_std = np.std(rho)
	
	if rho_std != 0:
		rho_hat = (rho-np.mean(rho))/np.std(rho)
	else:
		rho_hat = (rho-np.mean(rho))

	bic_hat = _standardise_bic(bic)
	VE = 0
	for itr_x in range(bic.shape[0]):
		for itr_y in range(bic.shape[1]):
			VE += abs(bic_hat[itr_x,itr_y] - rho_hat[itr_y])

	VE /= bic.shape[0]*bic.shape[1]
	return VE

def VEt(bic):
	"""
	source:

	Parameters:
		:param bic: numpy 2D array

	Output:
		:return : Virtual Error Transpose Score

	"""
	if bic.shape[0]<=1 or bic.shape[1]<=1:
		return 0.0
	return VE(np.transpose(bic))

def ASR(bic):
	"""
	source:

	Parameters:
		:param bic: numpy 2D array

	Output:
		:return : Average Spearman's Rho

	"""
	if bic.shape[0]<=1 or bic.shape[1]<=1:
		return 0.0

	spearman_rows = 0
	spearman_cols = 0
	for itr_x in range(bic.shape[0]-1):
		for itr_y in range(bic.shape[0]):
			spearman_rows += spearman(bic[itr_x,:],bic[itr_y,:])
	for itr_x in range(bic.shape[1]-1):
		for itr_y in range(bic.shape[1]):
			spearman_cols += spearman(bic[:,itr_x],bic[:,itr_y])

	spearman_rows /= bic.shape[0]*(bic.shape[0]-1)
	spearman_cols /= bic.shape[1]*(bic.shape[1]-1)
	value = 2*max(spearman_cols,spearman_rows)
	return value

def spearman(X,Y):
	"""
	Parameters:
		:param X: numpy vector
		:param Y: numpy vector of same size as X

	Output:
		:return : Spearman coefficient value 
	"""
	rankX = rankdata(X)
	rankY = rankdata(Y)
	m = len(X)
	coef = 6.0/(m**3 -m)
	value = 0
	for itr in range(m):
		value += (rankX[itr]-rankY[itr])**2
	return 1 - value*coef


def _standardise_bic(bic):
	"""
	Source:
	Pontes, B., Girldez, R., & Aguilar-Ruiz, J. S. (2015). Quality measures for gene expression biclusters. PloS one, 10(3), e0115497.

	Parameters:
		:param bic: numpy 2D array

	Output:
		:return  bic: Standardised bicluster by subtracting mean and dividing by standard deviation 
	"""
	bic_t = np.copy(bic)
	for itr_x in range(bic.shape[0]):
		g = bic[itr_x,:]
		std = np.std(g)
		if std !=0:
			bic_t[itr_x,:] = (g - np.mean(g))/std
		else:
			bic_t[itr_x,:] = (g - np.mean(g))

	return bic_t
