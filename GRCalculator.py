import sympy as sym

class metric:
    def __init__(self, matrix, variables):
        self.set_matrix(matrix)
        self.set_vars(variables)

    def set_matrix(self, matrix):
        self.__matrix = matrix

    def get_matrix(self):
        return self.__matrix

    def set_vars(self, variables):
        self.__variables = [x for x in variables]

    def get_var(self, i):
        return self.__variables[i]

    def get_elm(self, x, y):
        return self.__matrix[int(x*4 + y)]

    def get_inv(self):
        return self.__matrix.inv()

    def get_inv_elm(self, x, y):
        return self.__matrix.inv()[int(x*4 + y)]

class christoffel:
    def __init__(self, metric):
        self.set_metric(metric)

    def set_metric(self, metric):
        self.__metric = metric

    def get_metric(self):
        return self.__metric

    def solve(self, alpha, mu, nu):
        return (1/2)*(sum(self.__metric.get_inv_elm(alpha, i) * sym.diff(self.__metric.get_elm(nu, i), self.__metric.get_var(mu)) for i in range(4))
                    + sum(self.__metric.get_inv_elm(alpha, i) * sym.diff(self.__metric.get_elm(i, mu), self.__metric.get_var(nu)) for i in range(4))
                    + sum(self.__metric.get_inv_elm(alpha, i) * sym.diff(self.__metric.get_elm(mu, nu), self.__metric.get_var(i)) for i in range(4)) * -1
                )

class riemann:
    def __init__(self, christoffel):
        self.set_christoffel(christoffel)
        self.set_metric(christoffel.get_metric())

    def set_christoffel(self, christoffel):
        self.__christoffel = christoffel

    def get_christoffel(self):
        return self.__christoffel

    def set_metric(self, metric):
        self.__metric = metric
        self.__christoffel.set_metric(metric)

    def get_metric(self):
        return self.__metric

    def solve(self, alpha, beta, mu, nu):
        return (
            sym.diff(self.__christoffel.solve(alpha, beta, nu), self.__metric.get_var(mu)) - sym.diff(self.__christoffel.solve(alpha, beta, mu), self.__metric.get_var(nu))
            + sum(self.__christoffel.solve(alpha, mu, i)*self.__christoffel.solve(i, beta, nu) for i in range(4))
            - sum(self.__christoffel.solve(alpha, nu, i)*self.__christoffel.solve(i, beta, mu) for i in range(4))
        )

