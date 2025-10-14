from matplotlib import pyplot as plt


fig, ax = plt.subplots(1)

ax.plot(heat, T_hot, "o-", color="red")
ax.plot(heat, T_cold, "o-", color="blue")

ttd_u_location = (heat[1] * (1 + 1 / 15), T_cold[1]), (heat[1] * (1 + 1 / 15), T_hot[1])
bracket = mpatches.FancyArrowPatch(
    *ttd_u_location,
    connectionstyle=f"bar,angle=0,fraction=0",
    arrowstyle='-',                      # no arrowhead
    linewidth=2, color='gray'
)
ax.add_patch(bracket)

offset = heat[1] / 15     # absolute horizontal offset from the line
bar_len = heat[1] / 15    # absolute size of the bracket ends (in data units)
ax.plot([heat[1] + offset - bar_len/2, heat[1] + offset + bar_len/2], [T_hot[1], T_hot[1]], color='gray', lw=2)
ax.plot([heat[1] + offset - bar_len/2, heat[1] + offset + bar_len/2], [T_cold[1], T_cold[1]], color='gray', lw=2)

# Label
ax.text(heat[1] * 1.1, (T_hot[1] + T_cold[1]) / 2, r'$\Delta T = \text{ttd_u}$', va='center', color='gray', rotation=0)

ttd_l_location = (0 - offset, T_cold[0]), (0 - offset, T_hot[0])
bracket = mpatches.FancyArrowPatch(
    *ttd_l_location,
    connectionstyle=f"bar,angle=0,fraction=0",
    arrowstyle='-',                      # no arrowhead
    linewidth=2, color='gray'
)
ax.add_patch(bracket)

# Label
ax.text(heat[0] - offset * 5, (T_hot[0] + T_cold[0]) / 2, r'$\Delta T = \text{ttd_l}$', va='center', color='gray', rotation=0)
ax.plot([heat[0] - offset - bar_len/2, heat[0] - offset + bar_len/2], [T_hot[0], T_hot[0]], color='gray', lw=2)
ax.plot([heat[0] - offset - bar_len/2, heat[0] - offset + bar_len/2], [T_cold[0], T_cold[0]], color='gray', lw=2)
ax.set_xbound([- 6 * offset, heat[1] + 6 * offset])
