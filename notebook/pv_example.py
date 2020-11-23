from pvtrace import *
lsc = LSC((5.0, 5.0, 1.0))  # size in cm
lsc.show()                  # open visualiser
lsc.simulate(100)           # emit 100 rays
lsc.report()                # print report