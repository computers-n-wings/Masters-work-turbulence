import matplotlib.pyplot as plt
plt.rcParams["font.family"] = 'serif'
import math as m
import numpy as np
import matplotlib.font_manager as font_manager
import pylab as plot
params = {'legend.fontsize': 20,
          'legend.linewidth': 2}
plot.rcParams.update(params)


################################# Upload Data ###################################
#################################################################################

data1 = np.loadtxt("CHANNEL_DATA.txt")
y1 = data1[:,0]
u1 = data1[:,1]
v1 = data1[:,2]
w1 = data1[:,3]
uu1 = data1[:,4]
vv1 = data1[:,5]
ww1 = data1[:,6]

data2 = np.loadtxt("MOSER_KIM_MANSOUR.txt")
y2 = data2[:, 0]
yplus2 = data2[:, 1]
uuprime2 = data2[:, 2]
vvprime2 = data2[:, 3]
wwprime2 = data2[:, 4]
uvprime2 = data2[:, 5]
uwprime2 = data2[:, 6]
vwprime2 = data2[:, 7]

################################# Data Set 1 ####################################
#################################################################################

h = 1
Re_tau1 = 181.17

grad = np.gradient(u1, 0.005556)
Ub1 = np.trapz(u1, y1)/(2*h)

nu1 = (grad[0]*h**2)/Re_tau1**2
Rb1 = Ub1*h/nu1

u_tau = m.sqrt(nu1*grad[0])
delta_nu = nu1/u_tau

uplus1 = (1/u_tau)*u1
yplus1 = (1/delta_nu)*y1

u_ub = (1/Ub1)*u1

lny = np.zeros(64)
lob = np.zeros(64)

for i in range (0,64):
	lny[i] = m.log(yplus1[i+1])

[ka, B] = np.polyfit(lny[38:49], uplus1[39:50], 1)

for i in range (0,64):
	lob[i] = (lny[i] - B)/ka

print 1/ka

########################### Plots for Question 2 ###############################
################################################################################

fig = plt.figure(1)
fig.suptitle('Mean Flow Profile', fontsize=40, fontweight='bold')
label_size = 20
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.subplot(121)
plt.plot(yplus1[0:65], uplus1[0:65], 'r-o')
plt.xscale('log')
plt.xlabel( r'$y_+$', fontsize=28)
plt.ylabel(r'$\bar{u}_+$', fontsize=28)
axes = plt.gca()
axes.set_xlim([min(yplus1[0:65]),max(yplus1[0:65])])
plt.grid(True)

plt.subplot(122)
plt.plot(y1[0:65], u_ub[0:65], 'b-o')
plt.xscale('log')
plt.xlabel( r'$\frac{y}{h}$', fontsize=40)
plt.ylabel(r'$\frac{\bar{u}}{U_b}$', fontsize=40)
axes = plt.gca()
axes.set_xlim([min(y1[0:65]),max(y1[0:65])])
plt.grid(True)

plt.show()

################################ Question 3 ###################################
###############################################################################

uuprime1 = np.zeros_like(uuprime2)
vvprime1 = np.zeros_like(vvprime2)
wwprime1 = np.zeros_like(wwprime2)

for i in range(len(uuprime2)):
	uuprime1[i] = (uu1[i]-u1[i]**2)/(u_tau**2)
	vvprime1[i] = (vv1[i]-v1[i]**2)/(u_tau**2)
	wwprime1[i] = (ww1[i]-w1[i]**2)/(u_tau**2)

############################## uu stresses ####################################
###############################################################################

fig = plt.figure(1)
fig.suptitle('Mean Normal Stress in u direction for Data Set 1 and Data Set 2' 
	, fontsize=40, fontweight='bold')
label_size = 20
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size

plt.subplot(121)
plt.plot(yplus1[0:65], uuprime1, linewidth=2, 
	label = 'data set 1')
plt.xscale('log')
plt.hold(True)
plt.plot(yplus2, uuprime2,linewidth=2, linestyle='--',
	label = 'data set 2')
plt.xlabel( r'$y_+$', fontsize=32)
plt.ylabel(r'$\langle \frac{{{u}^\prime}^2}{u_\tau^2} \rangle$', fontsize=40)
axes = plt.gca()
axes.set_xlim([min(yplus1[1:65]),max(yplus1[1:65])])
plt.legend(loc="upper left")
plt.grid(True)

plt.subplot(122)
plt.plot(y1[0:65], uuprime1, linewidth=2, 
	label = 'data set 1')
plt.xscale('log')
plt.hold(True)
plt.plot(y2, uuprime2,linewidth=2, linestyle='--',
	label = 'data set 2')
plt.xscale('log')
plt.xlabel( r'$\frac{y}{h}$', fontsize=40)
plt.ylabel(r'$\langle \frac{{{u}^\prime}^2}{u_\tau^2} \rangle$', fontsize=40)
axes = plt.gca()
axes.set_xlim([min(y1[1:65]),max(y1[1:65])])
plt.legend(loc="upper left")
plt.grid(True)

plt.show()

########################## vv stresses #######################################
##############################################################################

fig = plt.figure(1)
fig.suptitle('Mean Normal Stress in v direction for Data Set 1 and Data Set 2'
	, fontsize=40, fontweight='bold')
label_size = 20
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size

plt.subplot(121)
plt.plot(yplus1[0:65], vvprime1, linewidth=2, 
	label = 'data set 1')
plt.xscale('log')
plt.hold(True)
plt.plot(yplus2, vvprime2,linewidth=2, linestyle='--', 
	label = 'data set 2')
plt.xlabel( r'$y_+$', fontsize=32)
plt.ylabel(r'$\langle \frac{{{v}^\prime}^2}{u_\tau^2} \rangle$', fontsize=40)
axes = plt.gca()
axes.set_xlim([min(yplus1[1:65]),max(yplus1[1:65])])
plt.legend(loc="upper left")

plt.grid(True)

plt.subplot(122)
plt.plot(y1[0:65], vvprime1, linewidth=2, 
	label = 'data set 1')
plt.xscale('log')
plt.hold(True)
plt.plot(y2, vvprime2,linewidth=2, linestyle='--', 
	label = 'data set 2')
plt.xscale('log')
plt.xlabel( r'$\frac{y}{h}$', fontsize=40)
plt.ylabel(r'$\langle \frac{{{v}^\prime}^2}{u_\tau^2} \rangle$', fontsize=40)
axes = plt.gca()
axes.set_xlim([min(y1[1:65]),max(y1[1:65])])
plt.legend(loc="upper left")

plt.grid(True)

plt.show()

################################ ww stresses ##################################
###############################################################################
fig = plt.figure(1)
fig.suptitle('Mean Normal Stress in w direction for Data Set 1 and Data Set 2',
 fontsize=40, fontweight='bold')
label_size = 20
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size

plt.subplot(121)
plt.plot(yplus1[0:65], wwprime1, linewidth=2, 
	label = 'data set 1')
plt.xscale('log')
plt.hold(True)
plt.plot(yplus2, wwprime2,linewidth=2, linestyle='--', 
	label = 'data set 2')
plt.xlabel( r'$y_+$', fontsize=32)
plt.ylabel(r'$\langle \frac{{{w}^\prime}^2}{u_\tau^2} \rangle$', fontsize=40)
axes = plt.gca()
axes.set_xlim([min(yplus1[1:65]),max(yplus1[1:65])])
plt.legend(loc="upper left")

plt.grid(True)

plt.subplot(122)
plt.plot(y1[0:65], wwprime1, linewidth=2, 
	label = 'data set 1')
plt.xscale('log')
plt.hold(True)
plt.plot(y2, wwprime2,linewidth=2, linestyle='--', 
	label = 'data set 2')
plt.xscale('log')
plt.xlabel( r'$\frac{y}{h}$', fontsize=40)
plt.ylabel(r'$\langle \frac{{{w}^\prime}^2}{u_\tau^2} \rangle$', fontsize=40)
axes = plt.gca()
axes.set_xlim([min(y1[1:65]),max(y1[1:65])])
plt.legend(loc="upper left")

plt.grid(True)

plt.show()


# ################################ Question 4 ###################################
# ###############################################################################
uvprime1 = np.zeros_like(uvprime2)

for i in range(len(uvprime2)):
	uvprime1[i] = (nu1*grad[i])/u_tau**2 - (1-y1[i])

plt.figure()
label_size = 20
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.semilogx(yplus1[0:65], uvprime1, linewidth=2, 
	label = r'$\langle u^\prime v^\prime \rangle$' ' for data set 1')
plt.title('Normalised Reynolds Stresses', 
	fontsize=40, fontweight='bold')
plt.xlabel( r'$y_+$', fontsize=32)
plt.ylabel(r'$\langle \frac{u^\prime v^\prime}{u_\tau^2} \rangle$', fontsize=40)
plt.hold(True)
plt.semilogx(yplus2, uvprime2,linewidth=2, linestyle='--', 
	label = r'$\langle u^\prime v^\prime \rangle$' ' for data set 2')
plt.legend()
axes = plt.gca()
axes.set_xlim([min(yplus1[1:65]),max(yplus1[1:65])])
plt.grid(True)
plt.show()

plt.figure()
label_size = 20
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.semilogx(y1[0:65], uvprime1, linewidth=2, 
	label = r'$\langle u^\prime v^\prime \rangle$' ' for data set 1')
plt.title('Normalised Reynolds Stresses', fontsize=40, fontweight='bold')
plt.xlabel( r'$\frac{y}{h}$', fontsize=32)
plt.ylabel(r'$\langle \frac{u^\prime v^\prime}{u_\tau^2} \rangle$', fontsize=40)
plt.hold(True)
plt.semilogx(y2, uvprime2,linewidth=2, linestyle='--', 
	label = r'$\langle u^\prime v^\prime \rangle$' ' for data set 2')
plt.legend()
axes = plt.gca()
axes.set_xlim([min(y1[1:65]),max(y1[1:65])])
plt.grid(True)
plt.show()