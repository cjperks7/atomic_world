'''

Estimates the fractional change in density from ICRF pumpout on WEST

cjperks
Apr 10, 2024

'''

# Modules
from atomic_world.run_TOFU import run_TOFU_WEST as rtfw
'''
# Shot #
shot =59989
#shot =60001

davg = {
    't1': [3.8, 4.3],
    't2': [5.4, 6.0]
    }
'''
shot = 59991
davg = {
    't1': [3.5,4.3],
    't2': [6.2,7.2]
}


rtfw.get_pumpout_est(
    device = 'WEST',
    shot = shot,
    sp = 'Ar',
    cs = 16,
    dt = 0.7,
    davg = davg
    )

ddata = {
    '59976': [1.133, 24.14, -4.46, 500, '48.9'],
    '59978': [3.666, 59.67, 3.647, 500, '48.9'],
    '59979': [7.355, 41.13, 2.168, 500, '48.9'],
    '59981': [16.82, -47.9, -4.51, 500, '48.9'],
    '59984': [25.87, 73.97, 41.35, 500, '48.9'],
    '59985': [32.93, 69.04, 75.76, 500, '48.9'],
    '59986': [37.12, 74.06, 69.12, 500, '48.9'],
    '59988': [28.03, 69.71, 14.80, 500, '53.0'],
    '59989': [28.22, 45.27, -18.3, 500, '53.0'],
    '59990_1': [32.59, 29.74, -4.74, 100, '53.0'],
    '59990_2': [31.23, 91.42, -17.9, 500, '53.0'],
    '59991': [35.31, 73.53, 56.27, 540, 'both'],
    '59992': [37.10, 65.40, 29.13, 500, '48.9'],
    '59993': [28.77, 66.18, 39.74, 500, '48.9'],
    '59995': [29.72, 28.24, -12.1, 500, '48.9'],
    '60001': [14.94, 31.01, 20.93, 500, '53.0']
}

plt.rcParams.update({'font.size': 20})
fig2, ax2 = plt.subplots()
msize = 15
lwid = 2

cnt_l = 0
num_l = 0
cnt_4 = 0
num_4 = 0
cnt_5 = 0
num_5 = 0

for shot in ddata.keys():
    if ddata[shot][-1] == '48.9':
        mark = "rv"
    elif ddata[shot][-1] == '53.0':
        mark = "b^"
    elif ddata[shot][-1] == 'both':
        mark = "g*"

    if ddata[shot][0] < 20:
        cnt_l += 1
        num_l += ddata[shot][2]
    else:
        if mark in ['rv', 'g*']:
            cnt_4 += 1
            num_4 += ddata[shot][2]
        else:
            cnt_5 += 1
            num_5 += ddata[shot][2]

    ax2.plot(
        ddata[shot][0],
        ddata[shot][2],
        mark,
        markersize = msize
        )

ax2.plot(
    [0,20],
    [num_l/cnt_l, num_l/cnt_l],
    'k--',
    linewidth = lwid +2
    )
ax2.plot(
    [20,40],
    [num_4/cnt_4, num_4/cnt_4],
    'r--',
    linewidth = lwid +2
    )
ax2.plot(
    [20,40],
    [num_5/cnt_5, num_5/cnt_5],
    'b--',
    linewidth = lwid +2
    )


ax2.set_ylabel('proxy density reduction [%]')
ax2.set_xlabel('H/H+D [%]')
ax2.grid('on')
ax2.set_ylim(-30, 100)
ax2.set_xlim(0, 50)









plt.rcParams.update({'font.size': 20})
fig1, ax1 = plt.subplots()
ax2 = ax1.twinx()
msize = 15
lwid = 2

cnt_l = 0
num_l = 0
cnt_4 = 0
num_4 = 0
cnt_5 = 0
num_5 = 0

for shot in ddata.keys():
    if ddata[shot][-1] == '48.9':
        mark = "v"
    elif ddata[shot][-1] == '53.0':
        mark = "^"
    elif ddata[shot][-1] == 'both':
        mark = "x"

    if ddata[shot][0] < 20:
        cnt_l += 1
        num_l += ddata[shot][2]
    else:
        if mark in ['x', 'v']:
            cnt_4 += 1
            num_4 += ddata[shot][2]
        else:
            cnt_5 += 1
            num_5 += ddata[shot][2]


    ax1.plot(
        [ddata[shot][0], ddata[shot][0]],
        [ddata[shot][1], ddata[shot][2]],
        'k--',
        linewidth = lwid
        )

    ax1.plot(
        ddata[shot][0],
        ddata[shot][1],
        color = 'blue',
        marker = mark,
        markersize= msize
        )
    ax2.plot(
        ddata[shot][0],
        ddata[shot][2],
        color = 'red',
        marker = mark,
        markersize = msize
        )

ax2.plot(
    [0,20],
    [num_l/cnt_l, num_l/cnt_l],
    'r--',
    linewidth = lwid +2
    )
ax2.plot(
    [20,40],
    [num_4/cnt_4, num_4/cnt_4],
    'r--',
    linewidth = lwid +2
    )
ax2.plot(
    [20,40],
    [num_5/cnt_5, num_5/cnt_5],
    'r--',
    linewidth = lwid +2
    )
ax2.plot(
    20, num_4/cnt_4,
    color = 'black',
    marker = 'v',
    markersize = msize 
    )
ax2.plot(
    40, num_4/cnt_4,
    color = 'black',
    marker = 'v',
    markersize = msize 
    )
ax2.plot(
    20, num_5/cnt_5,
    color = 'black',
    marker = '^',
    markersize = msize 
    )
ax2.plot(
    40, num_5/cnt_5,
    color = 'black',
    marker = '^',
    markersize = msize 
    )



ax1.set_ylabel('brightness reduction [%]', color='blue')
ax2.set_ylabel('proxy density reduction [%]', color = 'red')
ax1.set_xlabel('H/H+D [%]')
ax1.grid('on')
ax1.set_ylim(-50, 100)
ax2.set_ylim(-50, 100)
ax1.set_xlim(0, 60)