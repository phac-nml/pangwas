#!/usr/local/env python3

import os
import pangwas

def pytest_addoption(parser):
    parser.addoption("--validate", action="store_true", dest="validate", help="True if observed values should be validated.")
    parser.addoption("--generate", action="store_true", dest="generate", help="True if expected values should be generated.")
    parser.addoption("--unsafe",   action="store_false", dest="safe",    help="Allows execution of yaml variables.")
    parser.addoption("--tags", help='Rule tags to run.', nargs='+', default=[])
    #parser.addoption("--config", help='Path to the yaml config, relative to the tests directory.', default="test_pangwas.yaml")
    #parser.addoption("--profile", type=str, dest="profile", help="Profile to run. (default: test)", default="test")
    #parser.addoption("--profile-dir", type=str, dest="profile_dir", help="Profile directory(default: tests/profiles)", default=os.path.join("tests", "profiles"))
