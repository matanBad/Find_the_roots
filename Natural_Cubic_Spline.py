def solve_linear_system(matrix, vector):
    n = len(matrix)

    for i in range(n):
        max_row = i
        for k in range(i + 1, n):
            if abs(matrix[k][i]) > abs(matrix[max_row][i]):
                max_row = k

        matrix[i], matrix[max_row] = matrix[max_row], matrix[i]
        vector[i], vector[max_row] = vector[max_row], vector[i]

        if matrix[i][i] == 0:
            print("Cannot solve system: zero pivot.")
            return None

        for k in range(i + 1, n):
            factor = matrix[k][i] / matrix[i][i]

            for j in range(i, n):
                matrix[k][j] -= factor * matrix[i][j]

            vector[k] -= factor * vector[i]

    solution = [0.0] * n

    for i in range(n - 1, -1, -1):
        total = 0.0

        for j in range(i + 1, n):
            total += matrix[i][j] * solution[j]

        solution[i] = (vector[i] - total) / matrix[i][i]

    return solution


def cubic_spline_interpolation(x_points, y_points, x_target):
    if x_target in x_points:
        return y_points[x_points.index(x_target)]

    if x_target < min(x_points) or x_target > max(x_points):
        print("The target x is outside the table range.")
        return None

    n = len(x_points)

    h = []
    for i in range(n - 1):
        h.append(x_points[i + 1] - x_points[i])

    matrix = [[0.0 for _ in range(n)] for _ in range(n)]
    vector = [0.0 for _ in range(n)]

    # Natural spline conditions: second derivative at endpoints is 0
    matrix[0][0] = 1.0
    matrix[n - 1][n - 1] = 1.0

    for i in range(1, n - 1):
        matrix[i][i - 1] = h[i - 1]
        matrix[i][i] = 2 * (h[i - 1] + h[i])
        matrix[i][i + 1] = h[i]

        vector[i] = 6 * (
            ((y_points[i + 1] - y_points[i]) / h[i]) -
            ((y_points[i] - y_points[i - 1]) / h[i - 1])
        )

    second_derivatives = solve_linear_system(matrix, vector)

    if second_derivatives is None:
        return None

    interval_index = None
    for i in range(n - 1):
        if x_points[i] <= x_target <= x_points[i + 1]:
            interval_index = i
            break

    if interval_index is None:
        print("Could not find a suitable interval.")
        return None

    i = interval_index
    x_i = x_points[i]
    x_next = x_points[i + 1]
    y_i = y_points[i]
    y_next = y_points[i + 1]
    h_i = x_next - x_i

    m_i = second_derivatives[i]
    m_next = second_derivatives[i + 1]

    term1 = m_i * ((x_next - x_target) ** 3) / (6 * h_i)
    term2 = m_next * ((x_target - x_i) ** 3) / (6 * h_i)
    term3 = (y_i - (m_i * h_i ** 2) / 6) * ((x_next - x_target) / h_i)
    term4 = (y_next - (m_next * h_i ** 2) / 6) * ((x_target - x_i) / h_i)

    return term1 + term2 + term3 + term4


def main():
    x_points = [1.0, 2.0, 3.0, 4.0]
    y_points = [1.0, 4.0, 9.0, 16.0]

    x_target = 2.5

    print("Cubic Spline Interpolation Program")
    print("----------------------------------")

    print("Table points:")
    for i in range(len(x_points)):
        print(f"({x_points[i]}, {y_points[i]})")

    print(f"\nTarget x value: {x_target}")

    result = cubic_spline_interpolation(x_points, y_points, x_target)

    print("\nResult:")
    if result is not None:
        print(f"Cubic spline interpolation result: y ≈ {result}")


if __name__ == "__main__":
    main()