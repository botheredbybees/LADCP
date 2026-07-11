"""CastParams.from_ldeo_nc — reading cast parameters from LDEO_IX outputs."""
import netCDF4
import numpy as np

from ladcp.pipeline import CastParams


def _write_ref_nc(path, *, with_sadcp: bool, n_sadcp: int = 10):
    ds = netCDF4.Dataset(str(path), "w")
    ds.GEN_Magnetic_deviation_deg = -30.51
    ds.LADCP_dn_conf_single_ping_acc = 0.35
    ds.LADCP_dn_conf_number_pings = 10.0
    ds.uship = 0.009
    ds.vship = -0.0006
    ds.createDimension("lat", 1)
    v = ds.createVariable("lat", "f8", ("lat",))
    v[:] = -28.974
    if with_sadcp:
        ds.createDimension("z_sadcp", n_sadcp)
        for name, vals in (
            ("z_sadcp", np.linspace(35.0, 815.0, n_sadcp)),
            ("u_sadcp", np.linspace(-0.25, -0.05, n_sadcp)),
            ("v_sadcp", np.linspace(-0.09, 0.13, n_sadcp)),
            ("uerr_sadcp", np.full(n_sadcp, 0.03)),
        ):
            var = ds.createVariable(name, "f8", ("z_sadcp",))
            var[:] = vals
    ds.close()


def test_from_ldeo_nc_reads_attrs(tmp_path):
    p = tmp_path / "ref.nc"
    _write_ref_nc(p, with_sadcp=False)
    params = CastParams.from_ldeo_nc(p)
    assert params.lat_deg == -28.974
    assert params.drot_deg == -30.51
    assert abs(params.superens_std_min - 0.35 / np.sqrt(10.0)) < 1e-12
    assert params.u_ship == 0.009
    assert params.sadcp_z is None


def test_from_ldeo_nc_reads_embedded_sadcp(tmp_path):
    p = tmp_path / "ref_sadcp.nc"
    _write_ref_nc(p, with_sadcp=True)
    params = CastParams.from_ldeo_nc(p)
    assert params.sadcp_z is not None
    assert params.sadcp_z.shape == (10,)
    assert params.sadcp_u[0] == -0.25
    np.testing.assert_allclose(params.sadcp_err, 0.03)


def test_from_ldeo_nc_skips_all_nan_sadcp(tmp_path):
    p = tmp_path / "ref_nan.nc"
    ds = netCDF4.Dataset(str(p), "w")
    ds.GEN_Magnetic_deviation_deg = 0.0
    ds.createDimension("lat", 1)
    ds.createVariable("lat", "f8", ("lat",))[:] = 0.0
    ds.createDimension("z_sadcp", 5)
    for name in ("z_sadcp", "u_sadcp", "v_sadcp"):
        ds.createVariable(name, "f8", ("z_sadcp",))[:] = np.nan
    ds.close()
    params = CastParams.from_ldeo_nc(p)
    assert params.sadcp_z is None
