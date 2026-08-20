def test_package_is_importable() -> None:
    import hap_counter

    assert hap_counter.__doc__
