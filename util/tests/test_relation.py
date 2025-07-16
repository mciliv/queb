from unittest import TestCase
from collections import namedtuple
import itertools

from util.relation import Relation


class TestSets(TestCase):
    def test_sets(self):
        Params = namedtuple('Params', ['function_type', 'domain', 'range'])
        sets = itertools.product([Relation.one_to_one, Relation.bijection],
                                 [{0}, {1}, {0, 1}],
                                 [{'\x00'}, {'\x01'}, {'\x00', '\x01'}])
        sets = list(map(Params._make, list(sets)))
        actuals = [s.function_type(Relation(
            s.domain, s.range, lambda o, c: chr(o) == c)) for s in sets]
        results = [True, False, True,
                   False, True, True,
                   False, False, True,
                   
                   True, False, False,
                   False, True, False,
                   False, False, True]
        for i, s in enumerate(sets):
            with self.subTest(f'{s} and {results[i]}'):
                self.assertEqual(actuals[i], results[i])
