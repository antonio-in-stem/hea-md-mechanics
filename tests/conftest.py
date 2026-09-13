from pathlib import Path
import pytest
from hea_md.campaign import calculate

@pytest.fixture(scope='session')
def root():return Path(__file__).resolve().parents[1]

@pytest.fixture(scope='session')
def campaign_result(root):return calculate(root)
