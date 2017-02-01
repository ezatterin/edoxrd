


def plot_prof_compare(sample, sample2, d, d2, scale, L, H, win,polyorder):

        h,l,I = ut.read_rsm_data(sample, d, scale)
        h2, l2, I2 = ut.read_rsm_data(sample2, d2, scale)

        a = abs(l-L)
        b = abs(h-H)
        a2 = abs(l2-L)
        b2 = abs(h2-H)

        mask_a = ((a > 0) & (a < 1e-4))
        mask_b = ((b > 0) & (b < 1e-4))
        mask_a2 = ((a2 > 0) & (a2 < 1e-4))
        mask_b2 = ((b2 > 0) & (b2 < 1e-4))

        lst_a, lst_b = [], []
        lst_a2, lst_b2 = [], []

        for number in a[mask_a]:
            lst_a.append([(int(np.where(a==number)[0])),int(np.where(a==number)[1])])
        for number in b[mask_b]:
            lst_b.append([(int(np.where(b==number)[0])),int(np.where(b==number)[1])])
        for number in a2[mask_a2]:
            lst_a2.append([(int(np.where(a2==number)[0])),int(np.where(a2==number)[1])])
        for number in b2[mask_b2]:
            lst_b2.append([(int(np.where(b2==number)[0])),int(np.where(b2==number)[1])])

        idx_l = np.array(lst_a)
        idx_h = np.array(lst_b)
        idx_l2 = np.array(lst_a2)
        idx_h2 = np.array(lst_b2)

        I_h = I[idx_l[:,0], idx_l[:,1]][::-1]
        I_l = I[idx_h[:,0], idx_h[:,1]][::-1]
        I_h2 = I2[idx_l2[:,0], idx_l2[:,1]][::-1]
        I_l2 = I2[idx_h2[:,0], idx_h2[:,1]][::-1]

        h_h = h[idx_l[:,0], idx_l[:,1]][::-1]
        l_l = l[idx_h[:,0], idx_h[:,1]][::-1]
        h_h2 = h2[idx_l2[:,0], idx_l2[:,1]][::-1]
        l_l2 = l2[idx_h2[:,0], idx_h2[:,1]][::-1]

        # Plotting
        fig, ax = plt.subplots(1,4, figsize=(15,5))

        ax[0].pcolormesh(h,l,I); ax[0].set_ylabel('L'); ax[0].set_xlabel('H')
        ax[0].plot([h.min(),h.max()],[L,L],'r-')
        ax[0].plot([H, H], [l.min(),l.max()],'r-')
        ax[0].set_title(sample, fontsize=9)

        ax[1].pcolormesh(h2,l2,I2); ax[1].set_ylabel('L'); ax[1].set_xlabel('H')
        ax[1].plot([h2.min(),h2.max()],[L,L],'r-')
        ax[1].plot([H, H], [l2.min(),l2.max()],'r-')
        ax[1].set_title(sample2, fontsize=9)


        ax[2].plot(h_h, I_h, label=sample); ax[2].set_ylabel(r'$I$'); ax[2].set_xlabel('H')
        ax[2].plot(h_h2, I_h2, label=sample2); ax[2].set_ylabel(r'$I$'); ax[2].set_xlabel('H')
        plt.legend(fontsize=9)

        ax[3].plot(l_l, I_l, label=sample); ax[3].set_ylabel('$I$'); ax[3].set_xlabel('L')
        ax[3].plot(l_l2, I_l2, label=sample2); ax[3].set_ylabel('$I$'); ax[3].set_xlabel('L')
        plt.legend(fontsize=9)

        plt.tight_layout()
        plt.show()
