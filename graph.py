import matplotlib.pyplot as plt

# w = 0.8
errors1 = [5.48548, 1.12684, 0.270312, 0.0673144, 0.0171911, 0.00431299, 0.00106906, 0.000262786, 6.42253e-05, 1.62955e-05, 4.21592e-06, 1.0834e-06, 2.7705e-07, 7.05938e-08, 1.79402e-08, 4.55042e-09]
# w = 1
errors2 = [6.89173, 0.587316, 0.0216834, 0.00116604, 5.37742e-05, 2.37636e-06, 4.04162e-08, 1.18278e-09]

# w = 1.2
errors3 = [8.31165, 2.12583, 0.520333, 0.121019, 0.0313022, 0.0074822, 0.00180017, 0.000416651, 9.32689e-05, 2.02356e-05, 4.25313e-06, 8.63481e-07, 1.79608e-07, 3.93226e-08, 8.42911e-09]

n = len(errors1)

plt.xlabel("Iterations number")
plt.ylabel("Error")

plt.yscale('log')

plt.plot(range(len(errors1)), errors1)
plt.plot(range(len(errors2)), errors2)
plt.plot(range(len(errors3)), errors3)

plt.show()