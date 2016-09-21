import unittest
import os, sys, subprocess

class TestBasicUseCases(unittest.TestCase):

    def setUp(self):
        self.python_path = "/usr/local/bin/python"
        self.igv_plotter_path = os.path.dirname(__file__) + "/../bin/igv_plotter"

    def test_igv_plotter_exit_code(self):
        """Integration test for igv_plotter script"""
        exit_code = os.system(" ".join([self.python_path, self.igv_plotter_path, "1:12345"]))

        self.assertEqual(exit_code, 0)
