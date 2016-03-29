using CollectiveSpins

using PyCall
@pyimport matplotlib.pyplot as plt

function plot_geometry_xy(positions)
    for (i,p)=enumerate(positions)
        plt.text(p[1], p[2], "($i)", color="black", va="center", ha="center")
    end
    xmin = minimum([p[1] for p=positions])
    xmax = maximum([p[1] for p=positions])
    ymin = minimum([p[2] for p=positions])
    ymax = maximum([p[2] for p=positions])
    dx = xmax - xmin
    dy = ymax - ymin
    marginx = max(dx/5, 0.1)
    marginy = max(dy/5, 0.1)
    plt.xlim([xmin-marginx, xmax+marginx])
    plt.ylim([ymin-marginy, ymax+marginy])
end

d = 0.1
geomN = 5
positions = geometry.hexagonal(d; Nr=geomN)
positions = geometry.chain(d, 2*geomN+1)
plot_geometry_xy(positions)

plt.show()