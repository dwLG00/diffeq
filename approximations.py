import matplotlib.pyplot as plt
import math

def function_plot(f, start, end, n):
	x = [start + i*(end - start)/(n + 1) for i in range(n + 1)]
	l = [f(v) for v in x]
	return x, l

def riemann(f, start, subdiv):
	'''Generates integral function from f with a starting point and number of subdivisions'''

	'''

	def inner_function(x):
		if x == start:
			return 0
		elif x > start:
			n = (x - start)/incr
			s = 0
			for i in range(math.floor(n)):
				s += (f(start + (i+1)*incr) + f(start + i*incr)) * incr / 2
			s += (f(start + math.floor(n)*incr) + f(x)) * (n - math.floor(n)) / 2
			return s
		elif x < start:
			n = (start - x)/incr
			s = 0
			for i in range(math.floor(n)):
				s += -(f(start - (i+1)*incr) + f(start - i*incr)) * incr / 2
			s += -(f(start - math.floor(n)*incr) + f(x)) * (n - math.floor(n)) / 2
			return s

	return inner_function
	'''

	def inner_function(x):
		width = abs(x - start)
		incr = width / subdiv

		if x == start:
			return 0

		s = 0

		if x > start:
			for i in range(subdiv):
				cval = f(start + i*incr)
				nextval = f(start + (i+1)*incr)
				trapezoid = (cval + nextval) * incr / 2
				s += trapezoid
		elif x < start:
			for i in range(subdiv):
				cval = f(x + i*incr)
				nextval = f(x + (i+1)*incr)
				trapezoid = (cval + nextval) * incr / 2
				s += trapezoid
			s *= -1

		return s

	return inner_function

def flatten(f, a, b, slices):
	incr = (b - a)/slices
	l = [f(a + x*incr) for x in range(slices)]
	l.append(f(b)) #We store slices+1 values


	def func(x):
		idx = (x - a)/incr
		left = math.floor(idx)
		right = math.ceil(idx)
		decimal_portion = idx % 1
		return l[left] * (1-decimal_portion) + l[right] * decimal_portion

	return func

def picard(F, t0, x0, steps=5, riemann_subdiv=100, flatten_range=5):
	#F(t, x)

	flatten_slices = 2*riemann_subdiv

	xfuncs = [lambda t: x0]

	for i in range(steps):
		xfunc = xfuncs[-1]

		new_F = lambda t: F(t, xfunc(t))
		new_F_integral = riemann(new_F, t0, riemann_subdiv)
		xfunc_new = lambda t: x0 + new_F_integral(t)
		xfunc_new_flattened = flatten(xfunc_new, t0-flatten_range, t0+flatten_range, flatten_slices)
		xfuncs.append(xfunc_new_flattened)

	return xfuncs


def solve(alpha, beta, gamma, lambda_, f, p1, p2, fineness=10):
	'''(alpha + beta x + gamma x^2)dy/dx + lambda y = f(x) s.t. it passes through p1 and p2'''

	x1, y1 = p1
	x2, y2 = p2

	dydx = lambda x, y: (f(x) - (lambda_ * y))/(alpha + beta*x + gamma*x**2)

	solution_points = []
	first_point = p1 if x1 < x2 else p2
	second_point = p1 if x1 > x2 else p2
	dx = (second_point[0] - first_point[0])/(fineness - 1)

	solution_points.append(first_point)
	current_point = first_point

	for i in range(fineness - 1):
		cx, cy = current_point
		cdydx = dydx(cx, cy)
		next_point = (cx + dx, cy + cdydx*dx)
		solution_points.append(next_point)
		current_point = next_point

	return solution_points

def plot(points):
	x, y = [p[0] for p in points], [p[1] for p in points]
	plt.plot(x, y)

def simu():
	from pprint import pprint
	alpha = 1
	beta = 2
	gamma = 1
	lambda_ = 1
	f = lambda x: 1
	p1, p2 = (0, 5), (10, 0)

	solution = solve(alpha, beta, gamma, lambda_, f, p1, p2, fineness=100)
	return solution

def main():
	from pprint import pprint
	solution = simu()
	pprint(solution)

if __name__ == "__main__":
	main()
