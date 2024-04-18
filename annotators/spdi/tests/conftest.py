import pytest

@pytest.fixture
def mocked_wgs_reader(mocker):
    # This fixture will be used to mock the wgs_reader for each test in this class
    mocked_wgs_reader = mocker.Mock()
    return mocked_wgs_reader