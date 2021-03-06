from interval_matrix import DefaultIntervalMatrix, Interval, TwoDimmensionalTaskMatrix, TwoDimmensionalTaskVector, Matrix, IntervalVector
from interval_sole import IntervalSole
from nonlinear_interval_sole import NonlinearIntervalSole, JacobianInterval, Task, Test

import numpy as np
from matplotlib import pyplot as plt


def test():
    plt.plot([1, 3, 3, 1, 1], [1, 1, 3, 3, 1])
    plt.show()

    m = TwoDimmensionalTaskMatrix(4, 5, Interval(2, 3))
    print(m.at(0,0).interval_boundaries())
    print(m.at(0,1).interval_boundaries())
    print(m.at(1,0).interval_boundaries())
    print(m.at(1,1).interval_boundaries())

    mat = Matrix(2)
    mat[[0,0]] = 1
    mat[[0,1]] = 3
    mat[[1,0]] = 2
    mat[[1,1]] = 4

    im2 = mat.interval_matrix()

    print('Interval from dot')
    print(im2.at(0,0).interval_boundaries())
    print(im2.at(0,1).interval_boundaries())
    print(im2.at(1,0).interval_boundaries())
    print(im2.at(1,1).interval_boundaries())

    m1 = TwoDimmensionalTaskMatrix(2, 6, Interval(1, 2))
    im3 = im2.sub_interval_matrix(m1)

    print('sub')
    print(im3.at(0,0).interval_boundaries())
    print(im3.at(0,1).interval_boundaries())
    print(im3.at(1,0).interval_boundaries())
    print(im3.at(1,1).interval_boundaries())

    print('////////////////////////')
    im = m.mul_matrix(mat)

    print(im.at(0,0).interval_boundaries())
    print(im.at(0,1).interval_boundaries())
    print(im.at(1,0).interval_boundaries())
    print(im.at(1,1).interval_boundaries())

    iv = IntervalVector(2)
    iv.set(0, Interval(1, 2))
    iv.set(1, Interval(2, 3))
    
    print('////////////////////////')
    ivv = m.mul_vector(iv)

    print(ivv.at(0).interval_boundaries())
    print(ivv.at(1).interval_boundaries())

    print('////////////////////////')
    iv1 = mat.mul_interval_vector(iv)
    print(iv1.at(0).interval_boundaries())
    print(iv1.at(1).interval_boundaries())

    i1 = Interval(-3, 5)
    i2 = Interval(2, 6)
    print(Interval.intersection(i1, i2).interval_boundaries())

    id = Matrix.identity(2)
    print(id.matrix)


def nonlinear():
    pass

def main():
    # test()

    sole = IntervalSole(TwoDimmensionalTaskMatrix(1, 2, Interval(1, 2)), TwoDimmensionalTaskVector(Interval(2, 4)))
    sole.solve(1e-16)

    # NonlinearIntervalSole.create(np.array([
    #         Test.f1, Test.f2
    #     ]), JacobianInterval.create(np.array([
    #         np.array([Test.J11, Test.J12]),
    #         np.array([Test.J21, Test.J22])
    #     ])), IntervalVector.create(np.array([
    #     Interval(0, 1), Interval(0, 1)
    #     ]))).solve()

    sole = NonlinearIntervalSole.create(np.array([
            Task.f1, Task.f2
        ]), JacobianInterval.create(np.array([
            np.array([Task.J11, Task.J12]),
            np.array([Task.J21, Task.J22])
        ])), IntervalVector.create(np.array([
            Interval(0.25, 4), Interval(0.25, 4)
        ])))

    sole.solve(1e-16)

    return


if __name__ == '__main__':
    main()
