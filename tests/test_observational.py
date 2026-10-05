import unittest
import pandas as pd
from kg_audit.observational import landmark_records


class LandmarkTests(unittest.TestCase):
    def test_later_order_stays_in_landmark_comparator(self):
        frame=pd.DataFrame({'stay_id':[1,2,3],'t0':['2020-01-01']*3,'pre_t0_statin':[None]*3,
             'first_post_statin':['2020-01-02','2020-01-03 00:30','2020-01-05'],
             'deathtime':[None]*3,'dod':[None]*3,'dischtime':['2020-01-10']*3})
        out=landmark_records(frame)
        self.assertEqual(out.treat.tolist(),[1,0,0])
        self.assertTrue(out.alive_48h.eq(1).all())

    def test_day_only_death_straddling_landmark_is_flagged(self):
        frame=pd.DataFrame({'stay_id':[1],'t0':['2020-01-01 12:00'],'pre_t0_statin':[None],
             'first_post_statin':[None],'deathtime':[None],'dod':['2020-01-03'],'dischtime':['2020-01-10']})
        self.assertEqual(landmark_records(frame).ambiguous_landmark.iloc[0],1)

    def test_endpoint_horizon_is_28_days_from_entry_not_landmark(self):
        frame=pd.DataFrame({'stay_id':[1,2],'t0':['2020-01-01']*2,'pre_t0_statin':[None]*2,
             'first_post_statin':[None]*2,'deathtime':['2020-01-28','2020-01-30'],
             'dod':[None]*2,'dischtime':['2020-02-10']*2})
        self.assertEqual(landmark_records(frame).event28.tolist(),[1,0])
