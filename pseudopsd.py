from PAScual.PAScual import discretepals, MELTlikeROI, fitpar, distributeinsets
from PAScual.SpecFiles import ASCIIfileloader
from PAScual.pyTaoEldrup import TE_radius
import scipy as S


def fitpar_list_generator(*args):
    """generator for lists of fitpars corresponding to the given args.
    All args except one should be scalars. The one that is not a scalar
    should be an array of values for the given argument"""
    i = 0
    while S.isscalar(args[i]):
        i += 1
    ret = list(args)
    for e in args[i]:
        ret[i] = e
        yield [fitpar(e) for e in ret]


def pseudopsd_tau4():
    """
    pseudopsd is (min(chi2)/chi2)*ity4
    """
    file_loader = ASCIIfileloader(hdrlines=4)

    exp_fname = 'data/PTMSP_dry000.dat'
    ofname = 'data/psd_tau4.txt'

    left_of_max= 5
    stopdat = 840

    expdata = file_loader.expdata(exp_fname)
    roi = MELTlikeROI(expdata, left_of_max=left_of_max, stopdat=stopdat)

    tau4 = S.arange(500, 4000, 1000)
    chi2 = S.empty_like(tau4, dtype='float64')

    taulist_gen = fitpar_list_generator(125, 369, 1632, tau4, 6876)
    itylist = [fitpar(e) for e in (20.0, 47.3, 2.7, 5.3, 27.5)]


    for i,taulist in enumerate(taulist_gen):
        bg = fitpar(expdata[roi[-10:]].mean(), name='bg', free=True)
        fwhm = fitpar(260, name='fwhm', free=True)
        c0 = fitpar(50, name='c0', free=True)
        psperchannel = 49.5

        dpp = discretepals(name='tmp', expdata=expdata, roi=roi,
                           taulist=taulist, itylist=itylist, bg=bg, fwhm=fwhm,
                           c0=c0, psperchannel=psperchannel)

        pals_sets = distributeinsets([dpp])
        pset = pals_sets[0]
        pset.localmin()

        chi2[i] = dpp.recalculate_chi2(forcecalc=True) / dpp.dof
        dpp.showreport()



    ity4 = itylist[3].val
    pseudopsd = ity4 * chi2.min()/chi2
    pore_diameter = S.array([2 * TE_radius(tau/1000.) for tau in tau4])

    S.savetxt(ofname, S.stack((tau4, chi2, pore_diameter, pseudopsd)).T,
              header="tau4\tchi2\tpore_diameter(nm)\tpseudopsd",
              delimiter='\t', fmt='%8g')


    print "\n\n********* results written in '%s' ********\n\n" % ofname


if __name__ == '__main__':
    pseudopsd_tau4()