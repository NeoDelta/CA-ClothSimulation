#version 330 core
in vec3 Normal;
in vec3 Position;
out vec4 color;

uniform vec3 cameraPos;

void main()
{             
    vec3 I = normalize(Position - cameraPos);
    color = vec4(0.5, 0.0, 0.0, 1.0);
}