import numpy as np

# Matrix equation solver-------------------------------------------------------
def gaussSeidel(matrix_A, column_B, tolerance):
# -----------------------------------------------------------------------------
# 170620 There was division by zero because initial X vector for guess was 
# initialized as zero vector. By using numpy.ones instead of numpy.zeros, the
# issued was solved.

	"""	input)
		- matrix_A (square matrix)
		- column_B (column vector)
		  in [A]{x} = {B} equation
		- tolerance (convergence criterion) in %
	format)
		- matrix_A: np.matrix
		- column_B: np.array
	output)
		- X: solution column(vector)"""

	matrixSize = matrix_A.shape

	if matrixSize[0] != matrixSize[1]:
		print "warning) input should be square matrix!"
	if matrixSize[1] != len(column_B):
		print "warning) # of matrix's column not matched with # of vect\
				or's row!"
	if tolerance <= 0:
		print "warning) tolerance value should be positive!"

	else:
		X_length = len(column_B)
		X = np.ones(X_length, dtype=complex)
		X_new = X.copy() ## mere substitution makes two array share 
		running = True   ## properties. So, copy() method's required.
		iter = 0

		#Iterative loop
		while running:
			iter += 1

			for n in range(X_length):
				X_new[n] = column_B[n]

				i_range = range(X_length)
				del i_range[n]
				for i in i_range:
					X_new[n] -= X_new[i] * matrix_A[n, i]
				X_new[n] /= matrix_A[n, n]

			#Test the convergence
			if convergTest(X, X_new, iter, tolerance):
				break
			else:
				X = X_new.copy()

		return X_new

# The function that tests the convergence of the solution----------------------
def convergTest(Xold, Xnew, iter, tolerance):
# -----------------------------------------------------------------------------
	matrix_eps = (Xnew - Xold) / Xnew
        max_eps = np.amax(matrix_eps) * 100

	if max_eps < tolerance:
		print "{0}: The solution has converged!".format(iter)
		print "eps: {0} %".format(max_eps)
		return True
	else:
		print "{0}: The solution is being converged...".format(iter)
		print "eps: {0} %".format(max_eps)
		return False
