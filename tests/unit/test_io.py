from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from imz2anndata.io import detect_imzml_spectrum_mode, load_imzml_records


class ImzMLModeDetectionTests(unittest.TestCase):
    def _write_temp_imzml(self, content: str) -> Path:
        handle = tempfile.NamedTemporaryFile("w", suffix=".imzML", delete=False, encoding="utf-8")
        with handle:
            handle.write(content)
        return Path(handle.name)

    def test_detects_centroid_mode(self):
        path = self._write_temp_imzml('<cvParam accession="MS:1000127" cvRef="MS" name="centroid spectrum"/>')
        try:
            self.assertEqual(detect_imzml_spectrum_mode(path), "centroid")
        finally:
            path.unlink(missing_ok=True)

    def test_detects_profile_mode(self):
        path = self._write_temp_imzml('<cvParam accession="MS:1000128" cvRef="MS" name="profile spectrum"/>')
        try:
            self.assertEqual(detect_imzml_spectrum_mode(path), "profile")
        finally:
            path.unlink(missing_ok=True)

    def test_returns_none_when_mode_is_missing(self):
        path = self._write_temp_imzml("<mzML></mzML>")
        try:
            self.assertIsNone(detect_imzml_spectrum_mode(path))
        finally:
            path.unlink(missing_ok=True)

    def test_prefers_parser_spectrum_mode_when_available(self):
        from unittest import mock

        parser = mock.Mock()
        parser.spectrum_mode = "profile"
        parser.m = None

        with mock.patch("imz2anndata.io.ImzMLParser", return_value=parser) as parser_cls:
            path = Path("/tmp/test.imzML")
            with mock.patch.object(Path, "exists", return_value=True):
                self.assertEqual(detect_imzml_spectrum_mode(path), "profile")
        parser_cls.assert_called_once_with(str(path), ibd_file=None)

    def test_falls_back_to_text_search_when_parser_fails(self):
        from unittest import mock

        path = self._write_temp_imzml('<cvParam accession="MS:1000127" cvRef="MS" name="centroid spectrum"/>')
        try:
            with mock.patch("imz2anndata.io.ImzMLParser", side_effect=RuntimeError("boom")):
                self.assertEqual(detect_imzml_spectrum_mode(path), "centroid")
        finally:
            path.unlink(missing_ok=True)

    def test_load_imzml_records_returns_mode_and_records_from_parser(self):
        from unittest import mock

        parser = mock.Mock()
        parser.spectrum_mode = "centroid"
        parser.coordinates = [(1, 2, 1)]
        parser.getspectrum.return_value = ([100.0, 200.0], [10.0, 20.0])
        parser.m = None

        with mock.patch("imz2anndata.io.ImzMLParser", return_value=parser):
            mode, records = load_imzml_records(Path("/tmp/test.imzML"))

        self.assertEqual(mode, "centroid")
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].coordinate.x, 1)
        self.assertEqual(records[0].coordinate.y, 2)
        self.assertEqual(records[0].signal.mz.tolist(), [100.0, 200.0])


if __name__ == "__main__":
    unittest.main()
