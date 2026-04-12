import unittest
import tempfile
import sys
from unittest.mock import MagicMock

# Mock numpy to allow importing project files in the offline sandbox
sys.modules['numpy'] = MagicMock()

from readfile import read_files

class TestReadFile(unittest.TestCase):
    def test_read_files_missing_file(self):
        """Test that read_files raises FileNotFoundError when files are missing."""
        with tempfile.TemporaryDirectory() as temp_dir:
            with self.assertRaises(FileNotFoundError):
                read_files(temp_dir)

if __name__ == '__main__':
    unittest.main()