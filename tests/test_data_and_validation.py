"""Data, centering and held-out response regression checks."""
import sys
from pathlib import Path
import unittest
import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'experiment/code'))
from revision_data import build_panel,load_sample,rule_windows
from revision_empirical import centered_loss,curve_validation,fit_power,resample_counts


class DataAndValidationTests(unittest.TestCase):
    def test_missing_week_is_not_adjacent_observation(self):
        frame=pd.DataFrame({'state':['AA']*5,'week_end':pd.to_datetime(['2020-01-04','2020-01-11','2020-01-25','2020-02-01','2020-02-08']),
                            'C':[10.,20.,30.,40.,50.],'coverage':[100.]*5})
        panel=build_panel(frame)
        self.assertFalse(panel.consecutive.any())

    def test_missing_count_not_imputed_zero(self):
        frame=pd.DataFrame({'state':['AA']*6,'week_end':pd.date_range('2020-01-04',periods=6,freq='7D'),
                            'C':[10.,20.,np.nan,40.,50.,60.],'coverage':[100.]*6})
        p=build_panel(frame)
        self.assertTrue(p.loc[2,'any_missing'])
        self.assertTrue(np.isnan(p.loc[2,'C_t']))
        self.assertTrue(np.isnan(p.loc[3,'C_prev1']))

    def test_frozen_primary_sample(self):
        self.assertEqual(len(load_sample('RA')),1748)
        self.assertEqual(len(load_sample('G')),1756)
        self.assertFalse(load_sample().any_missing.any())
        self.assertEqual(sum((b-a).days//7+1 for a,b in rule_windows()[0].values()),42)

    def test_global_normalization_unchanged_after_season_center(self):
        f=load_sample();loss,valid,center=centered_loss(f,'season')
        means=f.groupby('season').x.transform('mean').to_numpy()
        np.testing.assert_allclose(loss,((f.x.to_numpy()-means)/f.x.mean())**2,atol=1e-14)

    def test_testfold_response_cannot_change_training_center(self):
        f=load_sample();first=sorted(f.season.unique())[0]
        original=curve_validation(f,'season')['folds'][0]
        changed=f.copy();changed.loc[changed.season==first,'x']*=1.8
        altered=curve_validation(changed,'season')['folds'][0]
        self.assertEqual(original['training_center'],altered['training_center'])
        self.assertEqual(original['affine_fit']['a'],altered['affine_fit']['a'])
        self.assertEqual(original['power_fit']['exponent'],altered['power_fit']['exponent'])
        self.assertNotEqual(original['score_affine'],altered['score_affine'])

    def test_power_known_curve(self):
        d=np.geomspace(5,5000,100);y=.4*d**(-.38)
        fit=fit_power(y,d)
        self.assertTrue(fit['converged'])
        self.assertAlmostEqual(fit['exponent'],.38,places=6)
        self.assertAlmostEqual(fit['coefficient'],.4,places=6)

    def test_blocks_preserve_total_weeks(self):
        f=load_sample();out=resample_counts({'RA':f},np.random.default_rng(123),3)
        weights,weeks=out['RA']
        self.assertEqual(weeks.sum(),42)
        self.assertTrue((weights>=0).all())


if __name__=='__main__':unittest.main()
