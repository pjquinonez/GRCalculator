import sympy as sym

class metric:
    def __init__(self, matrix, variables):
        self.set_matrix(matrix)
        self.set_variables(variables)

    def set_matrix(self, matrix):
        self.__matrix = matrix

    def get_matrix(self):
        return self.__matrix

    def set_variables(self, variables):
        self.__variables = [x for x in variables]

    def get_variable(self, i):
        return self.__variables[i]

    def get_element(self, x, y):
        return self.__matrix[int(x*4 + y)]

    def get_inverse(self):
        return self.__matrix.inv()

    def get_inverse_element(self, x, y):
        return self.__matrix.inv()[int(x*4 + y)]

class christoffel:
    def __init__(self, metric):
        self.set_metric(metric)

    def set_metric(self, metric):
        self.__metric = metric

    def solve(self, alpha, mu, nu):
        return (1/2)*(self.__metric.get_inverse_element(alpha, 0) * sym.diff(self.__metric.get_element(nu, 0), self.__metric.get_variable(mu))
                    + self.__metric.get_inverse_element(alpha, 1) * sym.diff(self.__metric.get_element(nu, 1), self.__metric.get_variable(mu))
                    + self.__metric.get_inverse_element(alpha, 2) * sym.diff(self.__metric.get_element(nu, 2), self.__metric.get_variable(mu))
                    + self.__metric.get_inverse_element(alpha, 3) * sym.diff(self.__metric.get_element(nu, 3), self.__metric.get_variable(mu))

                    + self.__metric.get_inverse_element(alpha, 0) * sym.diff(self.__metric.get_element(0, mu), self.__metric.get_variable(nu))
                    + self.__metric.get_inverse_element(alpha, 1) * sym.diff(self.__metric.get_element(1, mu), self.__metric.get_variable(nu))
                    + self.__metric.get_inverse_element(alpha, 2) * sym.diff(self.__metric.get_element(2, mu), self.__metric.get_variable(nu))
                    + self.__metric.get_inverse_element(alpha, 3) * sym.diff(self.__metric.get_element(3, mu), self.__metric.get_variable(nu))

                    - self.__metric.get_inverse_element(alpha, 0) * sym.diff(self.__metric.get_element(mu, nu), self.__metric.get_variable(0))
                    - self.__metric.get_inverse_element(alpha, 1) * sym.diff(self.__metric.get_element(mu, nu), self.__metric.get_variable(1))
                    - self.__metric.get_inverse_element(alpha, 2) * sym.diff(self.__metric.get_element(mu, nu), self.__metric.get_variable(2))
                    - self.__metric.get_inverse_element(alpha, 3) * sym.diff(self.__metric.get_element(mu, nu), self.__metric.get_variable(3))
                )