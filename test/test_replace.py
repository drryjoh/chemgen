f = open('test_read.h.in','r')
scalar = "double"
with f as file:
    # Read the file content
    content = file.read()
new_content  = eval(f'f"""{content}"""')
print(new_content)
