
from matplotlib.image import NonUniformImage
import matplotlib.pyplot as plt
import numpy as np
import time
j=6
grid = np.load('sq' + str(j) + '.npy')
grid = grid[0,:,:]
print np.shape(grid)

#fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, figsize=(6,10))

#ax1.imshow(grid, extent=[0,100,0,1])
#ax1.set_title('Default')
#Z=np.array(((1,2,3,4,5),(4,5,6,7,8),(7,8,9,10,11)))
#im = plt.imshow(Z, cmap='hot')
#plt.colorbar(im, orientation='horizontal')
plt.ion()
plt.imshow(grid[:,:], extent=[0,1,0,1], aspect='auto', vmin=0, vmax=1000)
plt.set_cmap('gray')
plt.colorbar()
plt.show()
plt.draw()
plt.savefig("f" + str(j).zfill(4)  + ".png",dpi=100)
time.sleep(10)

#ax2.imshow(grid, extent=[0,100,0,1], aspect='auto')
#plt.draw()
#ax2.set_title('Auto-scaled Aspect')

#ax3.imshow(grid, extent=[0,100,0,1], aspect=100)
#ax3.set_title('Manually Set Aspect')

