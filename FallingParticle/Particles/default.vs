#version 330 core
layout (location = 0) in vec3 position;

out vec3 Normal;
out vec3 Position;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

void main()
{
    gl_Position = projection * view * model * vec4(position, 1.0f);
    Position = vec3(model * vec4(position, 1.0f));
}  