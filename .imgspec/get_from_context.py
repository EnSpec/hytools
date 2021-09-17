import json
import sys

if __name__ == '__main__':
    with open("_context.json") as f:
        ctx = json.loads(f.read())
    field_name = sys.argv[1]
    value = ctx.get(field_name)
    if value is not None:
        print(json.dumps(value)) if type(value) == dict else print(value)
    else:
        value = None
        for params in ctx.get("job_specification").get("params"):
            if params.get("name") == field_name:
                value = params.get("value")
                print(json.dumps(value)) if type(value) == dict else print(value)
        if value is None:
            raise Exception(f"Field {field_name} doesn't exist in context file.")