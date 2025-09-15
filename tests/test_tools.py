"""Test FUSOR tools."""

from pydantic import BaseModel, ValidationError

from fusor.tools import get_error_message


class _TestModel(BaseModel):
    field1: int
    field2: str


def test_get_error_message():
    """Test that get_error_message works correctly"""
    # test single error message
    try:
        _TestModel(field1="not_an_int", field2="valid_str")
    except ValidationError as e:
        error_message = get_error_message(e)
        assert "should be a valid integer" in error_message

    # test multiple error messages in one ValidationError
    try:
        _TestModel(field1="not_an_int", field2=123)
    except ValidationError as e:
        error_message = get_error_message(e)
        assert "should be a valid integer" in error_message
        assert "should be a valid string" in error_message
