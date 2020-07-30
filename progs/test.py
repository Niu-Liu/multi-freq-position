fig, ax1 = plt.subplots()

ax1.set_xlabel("Bin size", fontsize=12)
ax1.set_ylabel("Correlation", color="r", fontsize=12)
ax1.plot(bin_array, x_corr["r_ka"], color="r", lw=1, ls="-", label="Pearson")
ax1.plot(bin_array, x_corr["rs_ka"], color="r", lw=1, ls="--", label="Spearman")
ax1.plot(bin_array, x_corr["tau_ka"], color="r", lw=1, ls="-.", label="Kendall")
ax1.hlines(0.5, 1, 51, lw=1, color="r", ls=":")
ax1.tick_params(axis="y", labelcolor="r")
ax1.set_xlim([0, 50])
ax1.legend(loc="upper left")

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

ax2.set_ylabel("$p$-value", color="b", fontsize=12)
ax2.plot(bin_array, x_corr["rp_ka"], color="b", lw=1, ls="-", label="Pearson")
ax2.plot(bin_array, x_corr["rsp_ka"], color="b", lw=1, ls="--", label="Spearman")
ax2.plot(bin_array, x_corr["taup_ka"], color="b", lw=1, ls="-.", label="Kendall")

# Critical value of p-value
ax2.hlines(0.001, 1, 51, lw=1, color="b", ls=":")
ax2.tick_params(axis="y", labelcolor="b")

plt.title("XKa-to-Gaia normalized offset vs. Gaia", fontsize=15)
fig.tight_layout()
