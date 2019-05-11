from domain import Domain
from billiards import *
from matplotlib import pyplot as plt
from matplotlib import animation


fourer_path = input("Enter a path to a Fourier coefficient file: ")
circle = Domain()
circle.import_fourier(fourer_path)
orbit_length = int(input("Enter the period desired: "))
#import pdb; pdb.set_trace()
bounce_theta = 0
bounce_angle = generate_orbit(circle, orbit_length)


fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                     xlim=(-1, 1), ylim=(-1, 1))
ax.set_xticks([])
ax.set_yticks([])
xdata, ydata = [], []
ln, = plt.plot([], [], 'b.-', markersize=0.5)

def init():
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    theta = 0
    while theta < 2 * pi:
        xdata.append(circle.polar(theta)[0])
        ydata.append(circle.polar(theta)[1])
        ln.set_data(xdata, ydata)
        theta += 0.01
    ln.set_data(xdata, ydata)
    return ln,

def update(frame):
    global bounce_theta, bounce_angle
    next_bounce = bounce(circle, bounce_theta, bounce_angle)
    bounce_theta = next_bounce[0]
    bounce_angle = next_bounce[1]
    xdata.append(circle.polar(bounce_theta)[0])
    ydata.append(circle.polar(bounce_theta)[1])
    ln.set_data(xdata, ydata)
    return ln,

anim = animation.FuncAnimation(fig, update, init_func=init,
                               frames=200, interval=100, blit=True)
plt.show()