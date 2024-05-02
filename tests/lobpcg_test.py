
import warnings
import numpy as np
from scipy.linalg import (inv, eigh, cho_factor, cho_solve, cholesky, LinAlgError)
from scipy.sparse.linalg import LinearOperator
from scipy.sparse import issparse

__all__ = ["lobpcg"]


def _report_nonhermitian(M, name):
    """
    Report if `M` is not a Hermitian matrix given its type.
    """
    from scipy.linalg import norm

    md = M - M.T.conj()
    nmd = norm(md, 1)
    tol = 10 * np.finfo(M.dtype).eps
    tol = max(tol, tol * norm(M, 1))
    if nmd > tol:
        warnings.warn(
              f"Matrix {name} of the type {M.dtype} is not Hermitian: "
              f"condition: {nmd} < {tol} fails.",
              UserWarning, stacklevel=4
         )

def _as2d(ar):
    """
    If the input array is 2D return it, if it is 1D, append a dimension,
    making it a column vector.
    """
    if ar.ndim == 2:
        return ar
    else:  # Assume 1!
        aux = np.array(ar, copy=False)
        aux.shape = (ar.shape[0], 1)
        return aux


def _makeMatMat(m):
    if m is None:
        return None
    elif callable(m):
        return lambda v: m(v)
    else:
        return lambda v: m @ v


def _matmul_inplace(x, y, verbosityLevel=0):
    """Perform 'np.matmul' in-place if possible.

    If some sufficient conditions for inplace matmul are met, do so.
    Otherwise try inplace update and fall back to overwrite if that fails.
    """
    if x.flags["CARRAY"] and x.shape[1] == y.shape[1] and x.dtype == y.dtype:
        # conditions where we can guarantee that inplace updates will work;
        # i.e. x is not a view/slice, x & y have compatible dtypes, and the
        # shape of the result of x @ y matches the shape of x.
        np.matmul(x, y, out=x)
    else:
        # ideally, we'd have an exhaustive list of conditions above when
        # inplace updates are possible; since we don't, we opportunistically
        # try if it works, and fall back to overwriting if necessary
        try:
            np.matmul(x, y, out=x)
        except Exception:
            if verbosityLevel:
                warnings.warn(
                    "Inplace update of x = x @ y failed, "
                    "x needs to be overwritten.",
                    UserWarning, stacklevel=3
                )
            x = x @ y
    return x


def _applyConstraints(blockVectorV,blockVectorY):
    """Changes blockVectorV in-place."""
    #YBV = blockVectorBY.T @ blockVectorV
    #tmp = cho_solve(factYBY, YBV)
    #blockVectorV -= blockVectorY @ tmp
    for i in range(blockVectorV.shape[1]):
        for j in range(blockVectorY.shape[1]):
            blockVectorV[:, i] -= np.dot(blockVectorY[:, j], blockVectorV[:, i]) / np.dot(blockVectorY[:, j], blockVectorY[:, j]) * blockVectorY[:, j]
        blockVectorV[:, i] /= np.sqrt(np.dot(blockVectorV[:, i], blockVectorV[:, i]))


def _b_orthonormalize(blockVectorV):
    """in-place B-orthonormalize the given block vector using Cholesky."""
    VBV = blockVectorV.T @ blockVectorV

    # VBV is a Cholesky factor from now on...
    #VBV = cholesky(VBV, overwrite_a=True)
    #VBV = inv(VBV, overwrite_a=True)
    #blockVectorV = _matmul_inplace(blockVectorV, VBV)

    """ for i in range(blockVectorV.shape[1]):
        for j in range(i):
            blockVectorV[:, i] -= np.dot(blockVectorV[:, j], blockVectorV[:, i]) / np.dot(blockVectorV[:, j], blockVectorV[:, j]) * blockVectorV[:, j]
        blockVectorV[:, i] /= np.sqrt(np.dot(blockVectorV[:, i], blockVectorV[:, i])) """
    for i in range(blockVectorV.shape[1]):
        V_i = blockVectorV[:, i]
        norm = np.sqrt(np.dot(V_i, V_i))
        V_i /= norm
        for j in range(i+1, blockVectorV.shape[1]):
            V_j = blockVectorV[:, j]
            proj = np.dot(V_i, V_j)
            V_j -= proj * V_i
            

    return blockVectorV, blockVectorV, VBV


def _get_indx(_lambda, num, largest):
    """Get `num` indices into `_lambda` depending on `largest` option."""
    ii = np.argsort(_lambda)
    if largest:
        ii = ii[:-num - 1:-1]
    else:
        ii = ii[:num]

    return ii


def _handle_gramA_gramB_verbosity(gramA, gramB, verbosityLevel):
    if verbosityLevel:
        _report_nonhermitian(gramA, "gramA")
        _report_nonhermitian(gramB, "gramB")


def lobpcg(
    A,
    X,
    B=None,
    M=None,
    Y=None,
    tol=None,
    maxiter=20,
    largest=True,
    verbosityLevel=0,
    retLambdaHistory=False,
    retResidualNormsHistory=False,
    restartControl=20,
):
    
    blockVectorX = X.copy()
    bestblockVectorX = blockVectorX
    blockVectorY = Y
    residualTolerance = tol

    bestIterationNumber = maxiter

    sizeY = 0
    if blockVectorY is not None:
        sizeY = blockVectorY.shape[1]

    n, sizeX = blockVectorX.shape

    # Data type of iterates, determined by X, must be inexact
    
    if (residualTolerance is None) or (residualTolerance <= 0.0):
        residualTolerance = np.sqrt(np.finfo(blockVectorX.dtype).eps) * n
    A = A.copy()
    A = _makeMatMat(A)
    B = _makeMatMat(B)
    M = _makeMatMat(M)

    # Apply constraints to X.
    if blockVectorY is not None:
        _applyConstraints(blockVectorX, blockVectorY)


    # B-orthonormalize X.
    blockVectorX, _, _ = _b_orthonormalize(blockVectorX)
    
    # Compute the initial Ritz vectors: solve the eigenproblem.
    blockVectorAX = A(blockVectorX)
    gramXAX = blockVectorX.T @ blockVectorAX

    _lambda, eigBlockVector = eigh(gramXAX, check_finite=False)
    ii = _get_indx(_lambda, sizeX, largest)
    _lambda = _lambda[ii]
    print("Eigenvectors Before Iterating: ", eigBlockVector)

    eigBlockVector = np.asarray(eigBlockVector[:, ii])
    for i in range(eigBlockVector.shape[1]):
        if eigBlockVector[0, i] < 0:
            eigBlockVector[:, i] *= -1
    
    print("gramXAX Before Iterating: ", gramXAX)
    
    blockVectorX = _matmul_inplace(blockVectorX, eigBlockVector)
    blockVectorAX = _matmul_inplace(blockVectorAX, eigBlockVector)


    ##
    # Active index set.
    activeMask = np.ones((sizeX,), dtype=bool)

    ##
    # Main iteration loop.

    blockVectorP = None  # set during iteration
    blockVectorAP = None

    smallestResidualNorm = np.abs(np.finfo(blockVectorX.dtype).max)

    iterationNumber = -1
    restart = True
    forcedRestart = False
    explicitGramFlag = False
    counter = 0
    while iterationNumber < maxiter:
        assert( maxiter > -1)
        iterationNumber += 1
        aux = blockVectorX * _lambda[np.newaxis, :]

        blockVectorR = blockVectorAX - aux
        aux = np.sum(blockVectorR * blockVectorR, 0)
        residualNorms = np.sqrt(np.abs(aux))
        print(f"Residual Norms: {residualNorms}")

        residualNorm = np.sum(np.abs(residualNorms)) / sizeX

        """ if residualNorm < smallestResidualNorm:
            smallestResidualNorm = residualNorm
            bestIterationNumber = iterationNumber
            bestblockVectorX = blockVectorX
        elif residualNorm > 2**restartControl * smallestResidualNorm:
            assert(False)
            forcedRestart = True
            blockVectorAX = A(blockVectorX)
 """
        ii = np.where(residualNorms > residualTolerance, True, False)
        activeMask = activeMask & ii
        currentBlockSize = activeMask.sum()
        print(f"Current Block Size: {currentBlockSize}")

        if currentBlockSize == 0:
            break

        activeBlockVectorR = _as2d(blockVectorR[:, :])

        if iterationNumber > 0:
            activeBlockVectorP = _as2d(blockVectorP[:, :])
            activeBlockVectorAP = _as2d(blockVectorAP[:, :])
        
        if M is not None:
            # Apply preconditioner T to the active residuals.
            activeBlockVectorR = M(activeBlockVectorR)

        ##
        # Apply constraints to the preconditioned residuals.
        if blockVectorY is not None:
            _applyConstraints(activeBlockVectorR, blockVectorY)

        ##
        # B-orthogonalize the preconditioned residuals to X.
        #activeBlockVectorR = activeBlockVectorR - (blockVectorX @ (blockVectorX.T @ activeBlockVectorR))
        _applyConstraints(activeBlockVectorR, blockVectorX)
        ##
        # B-orthonormalize the preconditioned residuals.
        aux = _b_orthonormalize(activeBlockVectorR)
        activeBlockVectorR, _, _ = aux


        activeBlockVectorAR = A(activeBlockVectorR)

        if iterationNumber > 0:
            aux = _b_orthonormalize(activeBlockVectorP)
            #_applyConstraints(activeBlockVectorP, np.concatenate((blockVectorX, activeBlockVectorR), axis=1))
            activeBlockVectorP, _, _ = aux
            # Function _b_orthonormalize returns None if Cholesky fails
            if activeBlockVectorP is not None:
                #activeBlockVectorAP = _matmul_inplace(activeBlockVectorAP, invR)
                restart = forcedRestart
            else:
                restart = True

        ##
        # Perform the Rayleigh Ritz Procedure:
        # Compute symmetric Gram matrices:

        if activeBlockVectorAR.dtype == "float32":
            myeps = 1
        else:
            myeps = np.sqrt(np.finfo(activeBlockVectorR.dtype).eps)

        if residualNorms.max() > myeps and not explicitGramFlag:
            explicitGramFlag = False
        else:
            # Once explicitGramFlag, forever explicitGramFlag.
            explicitGramFlag = True

        if not restart:
            activeBlockVectorAP = A(activeBlockVectorP)
         
            ##END OF STUFF WE MIGHT NOT NEED ################################
            subspace = np.concatenate((blockVectorX, activeBlockVectorR, activeBlockVectorP), axis=1)
            gramXAX = np.triu(np.dot(subspace.T, A(subspace)))
            #if iterationNumber==1:
                #print("Reference: ", gramA)
                #print("Reference - Mine: ", np.abs(gramA - gramA2))
                #print("RP ortho: ", np.dot(activeBlockVectorR.T, activeBlockVectorP))
                #print("XR ortho: ", np.dot(blockVectorX.T, activeBlockVectorR))
                #print("XP ortho: ", np.dot(blockVectorX.T, activeBlockVectorP))
            subspace, _,_ = _b_orthonormalize(subspace)
            print(f"Subspace: {subspace[0,:]}")
            #print(f"Subspace Shape: {subspace.shape}, AS shape: {A(subspace).shape}")
            gramXAX = np.dot(subspace.T, A(subspace))

            try:
                _lambda, eigBlockVector = eigh(gramXAX,
                                               check_finite=False)

            except LinAlgError as e:
                # raise ValueError("eigh failed in lobpcg iterations") from e
                if verbosityLevel:
                    warnings.warn(
                        f"eigh failed at iteration {iterationNumber} \n"
                        f"with error {e} causing a restart.\n",
                        UserWarning, stacklevel=2
                    )
                # try again after dropping the direction vectors P from RR
                restart = True

        if restart:
            counter += 1
            assert(counter < 2)
            subspace = np.concatenate((blockVectorX, activeBlockVectorR), axis=1)
            subspace, _,_ = _b_orthonormalize(subspace)
            print(f"Subspace: {subspace[0,:]}")
            gramXAX = np.dot(subspace.T, A(subspace))
            gramB = np.dot(subspace.T, subspace)

            try:
                _lambda, eigBlockVector = eigh(gramXAX,
                                               check_finite=False)
            except LinAlgError as e:
                # raise ValueError("eigh failed in lobpcg iterations") from e
                warnings.warn(
                    f"eigh failed at iteration {iterationNumber} with error\n"
                    f"{e}\n",
                    UserWarning, stacklevel=2
                )
                break

        ii = _get_indx(_lambda, sizeX, largest)
        _lambda = _lambda[ii]
        eigBlockVector = eigBlockVector[:, ii]
        #Invert the sign of the eigenvectors if the first element is negative
        for i in range(eigBlockVector.shape[1]):
            if eigBlockVector[0, i] < 0:
                eigBlockVector[:, i] *= -1
        print("Eigenvalues: ", _lambda)
        print("Eigenvectors of Iteration: ", eigBlockVector, iterationNumber)
        if retLambdaHistory:
            lambdaHistory[iterationNumber + 1, :] = _lambda

        # Compute Ritz vectors.
        print(f"Iteration Number, Restart: {iterationNumber}, {restart}")
        if not restart:
            eigBlockVectorX = eigBlockVector[:sizeX]
            eigBlockVectorR = eigBlockVector[sizeX:
                                                sizeX + currentBlockSize]
            eigBlockVectorP = eigBlockVector[sizeX + currentBlockSize:]

            pp = np.dot(activeBlockVectorR, eigBlockVectorR)
            pp += np.dot(activeBlockVectorP, eigBlockVectorP)

            app = np.dot(activeBlockVectorAR, eigBlockVectorR)
            app += np.dot(activeBlockVectorAP, eigBlockVectorP)
            #print(f"BlockVectorX Prior to Scaling: {blockVectorX}")
            #print(f"BlockVectorR Used in Scaling: {activeBlockVectorR}")
            #print(f"BlockVectorP Used in Scaling: {activeBlockVectorP}")
        else:
            eigBlockVectorX = eigBlockVector[:sizeX]
            eigBlockVectorR = eigBlockVector[sizeX:]
            #print(f"BlockVectorR Used in Scaling: {activeBlockVectorR}")
            pp = np.dot(activeBlockVectorR, eigBlockVectorR)
            app = np.dot(activeBlockVectorAR, eigBlockVectorR)
        

        blockVectorX = np.dot(blockVectorX, eigBlockVectorX) + pp
        blockVectorAX = np.dot(blockVectorAX, eigBlockVectorX) + app

        blockVectorP, blockVectorAP = pp, app
        

        #print(f"BlockR: {activeBlockVectorR}")
        #print(f"BlockP: {blockVectorP}")
        #print(f"BlockAX: {blockVectorAX}")
    finalEigVects = eigBlockVector.copy()
    finalgramXAX = gramXAX.copy()
    if maxiter == -1:
        subspace = blockVectorX.copy()
    elif maxiter == 0:
        subspace = np.concatenate((blockVectorX.copy(), activeBlockVectorR.copy()), axis=1)
    else:
        subspace = np.concatenate((blockVectorX.copy(), activeBlockVectorR.copy(), blockVectorP.copy()), axis=1)

    aux = blockVectorX * _lambda[np.newaxis, :]

    blockVectorR = blockVectorAX - aux

    aux = np.sum(blockVectorR.conj() * blockVectorR, 0)
    residualNorms = np.sqrt(np.abs(aux))
    # Use old lambda in case of early loop exit.
    residualNorm = np.sum(np.abs(residualNorms)) / sizeX
    if residualNorm < smallestResidualNorm:
        smallestResidualNorm = residualNorm
        bestIterationNumber = iterationNumber + 1
        bestblockVectorX = blockVectorX

    if np.max(np.abs(residualNorms)) > residualTolerance:
        warnings.warn(
            f"Exited at iteration {iterationNumber} with accuracies \n"
            f"{residualNorms}\n"
            f"not reaching the requested tolerance {residualTolerance}.\n"
            f"Use iteration {bestIterationNumber} instead with accuracy \n"
            f"{smallestResidualNorm}.\n",
            UserWarning, stacklevel=2
        )

    if verbosityLevel:
        print(f"Final iterative eigenvalue(s):\n{_lambda}")
        print(f"Final iterative residual norm(s):\n{residualNorms}")

    blockVectorX = bestblockVectorX
    # Making eigenvectors "exactly" satisfy the blockVectorY constrains
    if blockVectorY is not None:
        _applyConstraints(blockVectorX, blockVectorY)

    # Making eigenvectors "exactly" othonormalized by final "exact" RR
    blockVectorAX = A(blockVectorX)
    if blockVectorAX.shape != blockVectorX.shape:
        raise ValueError(
            f"The shape {blockVectorX.shape} "
            f"of the postprocessing iterate not preserved\n"
            f"and changed to {blockVectorAX.shape} "
            f"after multiplying by the primary matrix.\n"
        )
    gramXAX = np.dot(blockVectorX.T, blockVectorAX)


    gramXBX = np.dot(blockVectorX.T, blockVectorX)

    gramXAX = (gramXAX + gramXAX.T) / 2
    gramXBX = (gramXBX + gramXBX.T) / 2
    try:
        _lambda, eigBlockVector = eigh(gramXAX,
                                       gramXBX,
                                       check_finite=False)
    except LinAlgError as e:
        raise ValueError("eigh has failed in lobpcg postprocessing") from e

    ii = _get_indx(_lambda, sizeX, largest)
    _lambda = _lambda[ii]
    eigBlockVector = np.asarray(eigBlockVector[:, ii])

    blockVectorX = np.dot(blockVectorX, eigBlockVector)
    blockVectorAX = np.dot(blockVectorAX, eigBlockVector)

    aux = blockVectorX * _lambda[np.newaxis, :]
    
    blockVectorR = blockVectorAX - aux

    aux = np.sum(blockVectorR.conj() * blockVectorR, 0)
    residualNorms = np.sqrt(np.abs(aux))

    if retLambdaHistory:
        lambdaHistory[bestIterationNumber + 1, :] = _lambda
    if retResidualNormsHistory:
        residualNormsHistory[bestIterationNumber + 1, :] = residualNorms

    if retLambdaHistory:
        lambdaHistory = lambdaHistory[
            : bestIterationNumber + 2, :]
    if retResidualNormsHistory:
        residualNormsHistory = residualNormsHistory[
            : bestIterationNumber + 2, :]

    if np.max(np.abs(residualNorms)) > residualTolerance:
        pass

    if verbosityLevel:
        print(f"Final postprocessing eigenvalue(s):\n{_lambda}")
        print(f"Final residual norm(s):\n{residualNorms}")

    if retLambdaHistory:
        lambdaHistory = np.vsplit(lambdaHistory, np.shape(lambdaHistory)[0])
        lambdaHistory = [np.squeeze(i) for i in lambdaHistory]
    if retResidualNormsHistory:
        residualNormsHistory = np.vsplit(residualNormsHistory,
                                         np.shape(residualNormsHistory)[0])
        residualNormsHistory = [np.squeeze(i) for i in residualNormsHistory]
    
    if retLambdaHistory:
        if retResidualNormsHistory:
            return _lambda, subspace, finalEigVects, finalgramXAX, lambdaHistory, residualNormsHistory
        else:
            return _lambda, subspace, finalEigVects, finalgramXAX, lambdaHistory
    else:
        if retResidualNormsHistory:
            return _lambda, subspace, finalEigVects, finalgramXAX, residualNormsHistory
        else:
            return _lambda, subspace, finalEigVects, finalgramXAX
