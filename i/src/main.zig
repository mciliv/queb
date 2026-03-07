const std = @import("std");
const ai = @import("ai.zig");

pub fn main() !void {
    const allocator = std.heap.page_allocator;

    // EXAMPLE 1: OpenAI
    const openai_key = std.process.getEnvVarOwned(allocator, "OPENAI_API_KEY") catch "";
    defer if (openai_key.len > 0) allocator.free(openai_key);
    if (openai_key.len == 0) {
        std.debug.print("Set OPENAI_API_KEY for remote OpenAI usage\n", .{});
    }

    var client = ai.Client.init(allocator, .{ .openai = .{
        .base_url = "https://api.openai.com",
        .api_key = openai_key,
    } });

    const resp1 = try client.chat(.{
        .model = "gpt-4.1",
        .messages = &[_]ai.Message{
            .{ .role = .user, .content = "Zig is cool, explain why." },
        },
    });
    std.debug.print("OpenAI: {s}\n", .{resp1.text});

    // EXAMPLE 2: Local Ollama
    var local_client = ai.Client.init(allocator, .{ .local = .{
        .base_url = "http://localhost:11434",
    } });
    const resp2 = try local_client.chat(.{
        .model = "llama3.2",
        .messages = &[_]ai.Message{
            .{ .role = .user, .content = "Hello from local" },
        },
    });
    std.debug.print("Local: {s}\n", .{resp2.text});
}
