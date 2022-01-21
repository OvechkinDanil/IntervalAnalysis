import numpy as np
import kaucherpy as kp
import matplotlib.pyplot as plt
import random

def get_Kaucher_mat(A_inf, A_sup):
  n = A_inf.shape[0]
  return np.array([[kp.Kaucher(A_inf[i][j], A_sup[i][j]) for j in range(n)] for i in range(n)])

def get_inv_Kaucher_mat(A):
  n = A.shape[0]
  A_inf = np.array([[A[i][j].lower for j in range(n)] for i in range(n)])
  A_sup = np.array([[A[i][j].upper for j in range(n)] for i in range(n)])
  return (A_inf, A_sup)

def get_Kaucher_vec(b_inf, b_sup):
  n = b_inf.shape[0]
  return np.array([kp.Kaucher(b_inf[i], b_sup[i]) for i in range(n)])

def get_inv_Kaucher_vec(b):
  n = b.shape[0]
  b_inf = np.array([b[i].lower for i in range(n)])
  b_sup = np.array([b[i].upper for i in range(n)])
  return (b_inf, b_sup)

def get_mid_Kaucher_vec(b):
  b_mid = np.zeros(b.shape[0])
  for i in range(b.shape[0]):
    b_mid[i] = (b[i].upper + b[i].lower) / 2
  return b_mid

def _get_sti_vec(vec_inf, vec_sup):
    return np.append(-vec_inf, vec_sup)


def _get_inv_sti_vec(vec):
    n = vec.shape[0] // 2
    vec_inf = -vec[:n]
    vec_sup = vec[n:]
    return (vec_inf, vec_sup)


def _get_sign_block_matrix(mat):
    mat_pos = mat.copy()
    mat_neg = mat.copy()
    mat_pos[mat_pos < 0] = 0
    mat_neg[mat_neg > 0] = 0
    return np.block([[mat_pos, -mat_neg],
                     [-mat_neg, mat_pos]])


def _get_x0(C_inf, C_sup, d_inf, d_sup):
    C_mid_blocked = _get_sign_block_matrix((C_sup + C_inf) / 2)
    d_sti = _get_sti_vec(d_inf, d_sup)
    return np.linalg.solve(C_mid_blocked, d_sti)


def _Gxk(C, y, d):
    y_inv = get_Kaucher_vec(*_get_inv_sti_vec(y))
    dot_prod = C.dot(y_inv)
    sti_prod = _get_sti_vec(*get_inv_Kaucher_vec(dot_prod))
    return sti_prod - _get_sti_vec(*get_inv_Kaucher_vec(d))


def _Dxk(A, x):
    n = A.shape[0]
    D = np.zeros((2 * n, 2 * n))

    for i in range(n):
        for j in range(n):
            xinf = x[j].lower
            xsup = x[j].upper
            ainf = A[i][j].lower
            asup = A[i][j].upper

            k = 0
            m = 0

            if ainf * asup > 0:
                if ainf > 0:
                    k = 0
                else:
                    k = 2
            else:
                if ainf < asup:
                    k = 1
                else:
                    k = 3

            if xinf * xsup > 0:
                if xinf > 0:
                    m = 1
                else:
                    m = 3
            else:
                if xinf <= xsup:
                    m = 2
                else:
                    m = 4

            case = 4 * k + m

            if case == 1:
                D[i][j] = ainf
                D[i + n][j + n] = asup
            elif case == 2:
                D[i][j] = asup
                D[i + n][j + n] = asup
            elif case == 3:
                D[i][j] = asup
                D[i + n][j + n] = ainf
            elif case == 4:
                D[i][j] = ainf
                D[i + n][j + n] = ainf
            elif case == 5:
                D[i][j + n] = ainf
                D[i + n][j + n] = asup
            elif case == 6:
                if ainf * xsup < asup * xinf:
                    D[i][j + n] = ainf
                else:
                    D[i][j] = asup
                if ainf * xinf > asup * xsup:
                    D[i + n][j] = ainf
                else:
                    D[i + n][j + n] = asup
            elif case == 7:
                D[i][j] = asup
                D[i + n][j] = ainf
            elif case == 8 or case == 9:
                D[i][j + n] = ainf
                D[i + n][j] = asup
            elif case == 10:
                D[i][j + n] = ainf
                D[i + n][j] = ainf
            elif case == 11:
                D[i][j + n] = asup
                D[i + n][j] = ainf
            elif case == 12:
                D[i][j + n] = asup
                D[i + n][j] = asup
            elif case == 13:
                D[i][j] = ainf
                D[i + n][j] = asup
            elif case == 14 or case == 15:
                D[i][j + n] = asup
                D[i + n][j + n] = ainf
            elif case == 16:
                if ainf * xinf > asup * xsup:
                    D[i][j] = ainf
                else:
                    D[i][j + n] = -asup
                if ainf * xsup < asup * xinf:
                    D[i + n][j + n] = ainf
                else:
                    D[i + n][j] = asup
    return D


def _get_sti_vec(vec_inf, vec_sup):
    return np.append(-vec_inf, vec_sup)


def _get_inv_sti_vec(vec):
    n = vec.shape[0] // 2
    vec_inf = -vec[:n]
    vec_sup = vec[n:]
    return (vec_inf, vec_sup)


def _get_sign_block_matrix(mat):
    mat_pos = mat.copy()
    mat_neg = mat.copy()
    mat_pos[mat_pos < 0] = 0
    mat_neg[mat_neg > 0] = 0
    return np.block([[mat_pos, -mat_neg],
                     [-mat_neg, mat_pos]])


def _get_x0(C_inf, C_sup, d_inf, d_sup):
    C_mid_blocked = _get_sign_block_matrix((C_sup + C_inf) / 2)
    d_sti = _get_sti_vec(d_inf, d_sup)
    return np.linalg.solve(C_mid_blocked, d_sti)


def _Gxk(C, y, d):
    y_inv = get_Kaucher_vec(*_get_inv_sti_vec(y))
    dot_prod = C.dot(y_inv)
    sti_prod = _get_sti_vec(*get_inv_Kaucher_vec(dot_prod))
    return sti_prod - _get_sti_vec(*get_inv_Kaucher_vec(d))


def _Dxk(A, x):
    n = A.shape[0]
    D = np.zeros((2 * n, 2 * n))

    for i in range(n):
        for j in range(n):
            xinf = x[j].lower
            xsup = x[j].upper
            ainf = A[i][j].lower
            asup = A[i][j].upper

            k = 0
            m = 0

            if ainf * asup > 0:
                if ainf > 0:
                    k = 0
                else:
                    k = 2
            else:
                if ainf < asup:
                    k = 1
                else:
                    k = 3

            if xinf * xsup > 0:
                if xinf > 0:
                    m = 1
                else:
                    m = 3
            else:
                if xinf <= xsup:
                    m = 2
                else:
                    m = 4

            case = 4 * k + m

            if case == 1:
                D[i][j] = ainf
                D[i + n][j + n] = asup
            elif case == 2:
                D[i][j] = asup
                D[i + n][j + n] = asup
            elif case == 3:
                D[i][j] = asup
                D[i + n][j + n] = ainf
            elif case == 4:
                D[i][j] = ainf
                D[i + n][j + n] = ainf
            elif case == 5:
                D[i][j + n] = ainf
                D[i + n][j + n] = asup
            elif case == 6:
                if ainf * xsup < asup * xinf:
                    D[i][j + n] = ainf
                else:
                    D[i][j] = asup
                if ainf * xinf > asup * xsup:
                    D[i + n][j] = ainf
                else:
                    D[i + n][j + n] = asup
            elif case == 7:
                D[i][j] = asup
                D[i + n][j] = ainf
            elif case == 8 or case == 9:
                D[i][j + n] = ainf
                D[i + n][j] = asup
            elif case == 10:
                D[i][j + n] = ainf
                D[i + n][j] = ainf
            elif case == 11:
                D[i][j + n] = asup
                D[i + n][j] = ainf
            elif case == 12:
                D[i][j + n] = asup
                D[i + n][j] = asup
            elif case == 13:
                D[i][j] = ainf
                D[i + n][j] = asup
            elif case == 14 or case == 15:
                D[i][j + n] = asup
                D[i + n][j + n] = ainf
            elif case == 16:
                if ainf * xinf > asup * xsup:
                    D[i][j] = ainf
                else:
                    D[i][j + n] = -asup
                if ainf * xsup < asup * xinf:
                    D[i + n][j + n] = ainf
                else:
                    D[i + n][j] = asup
    return D


def printMatrix(A):
    print("MATRIX")
    matr = np.array(A)
    for line in matr:
        index = 0
        line = np.array(line)
        result = ""
        for elem in line:
            index += 1
            if index != 18:
                result += f"{np.round(elem, 3)} & "

        result += " \\\ "
        print(result)


def subDiff2(A_inf, A_sup, b_inf, b_sup, relax=1, max_iter=150, eps=1e-5):
    C = get_Kaucher_mat(A_inf, A_sup)
    d = get_Kaucher_vec(b_inf, b_sup)

    # printMatrix(C)
    print("ITER = 0")
    x_prev = _get_x0(A_inf, A_sup, b_inf, b_sup)
    print(f"X(0) = {x_prev}")
    G = _Gxk(C, x_prev, d)
    D = _Dxk(C, get_Kaucher_vec(*_get_inv_sti_vec(x_prev)))
    print(f"Субградиент отображения: {np.linalg.solve(D, G)}")
    x_next = x_prev - np.linalg.solve(D, G)
    iter = 1

    print(f"X(ITER) = {x_next}")

    norms = []
    norms.append(np.linalg.norm(x_next - x_prev))
    print(f"isEscape? {np.linalg.norm(x_next - x_prev) >= eps}")
    print(f"|| x_next - x_prev || = {np.linalg.norm(x_next - x_prev)}")
    while iter <= max_iter and np.linalg.norm(x_next - x_prev) >= eps:
        x_prev = x_next
        G = _Gxk(C, x_prev, d)
        D = _Dxk(C, get_Kaucher_vec(*_get_inv_sti_vec(x_prev)))
        x_next = x_prev - np.linalg.solve(D, G)
        print(f"Субградиент отображения: {np.linalg.solve(D, G)}")
        print(f"X({iter + 1} = {x_next}")
        iter = iter + 1
        norms.append(np.linalg.norm(x_next - x_prev))

    return (_get_inv_sti_vec(x_next), iter, norms)
