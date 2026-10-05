import unittest
from kg_audit.splits import validation_partition


class SplitTests(unittest.TestCase):
    def test_edge_partition_is_disjoint_and_complete(self):
        pairs=[(f'd{i}',f'x{j}') for i in range(8) for j in range(3)]
        train,val=validation_partition(pairs,'random',7)
        self.assertFalse(set(train)&set(val));self.assertEqual(set(train)|set(val),set(pairs))

    def test_cold_validation_uses_unlabelled_drugs(self):
        pairs=[(f'd{i}',f'x{j}') for i in range(8) for j in range(3)]
        train,val=validation_partition(pairs,'compound_disjoint',7)
        self.assertFalse({c for c,d in train}&{c for c,d in val})
        self.assertEqual(set(train)|set(val),set(pairs))

    def test_small_unidentifiable_task_rejected(self):
        with self.assertRaises(ValueError):validation_partition([('d','x'),('d','y')],'compound_disjoint',7)
