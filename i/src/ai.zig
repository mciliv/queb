const std = @import("std");

pub const Role = enum { system, user, assistant };

pub const Message = struct {
    role: Role,
    content: []const u8,
};

pub const ChatRequest = struct {
    model: []const u8,
    messages: []const Message,
    temperature: ?f32 = null,
};

pub const ChatResponse = struct {
    text: []const u8,
};

pub const Backend = union(enum) {
    openai: OpenAI,
    local: Local, // local like Ollama
};

pub const OpenAI = struct {
    base_url: []const u8,
    api_key: []const u8,
};

pub const Local = struct {
    base_url: []const u8,
};

pub const Client = struct {
    allocator: std.mem.Allocator,
    backend: Backend,

    pub fn init(allocator: std.mem.Allocator, backend: Backend) Client {
        return Client{ .allocator = allocator, .backend = backend };
    }

    pub fn chat(self: *Client, req: ChatRequest) !ChatResponse {
        return switch (self.backend) {
            .openai => |cfg| openaiChat(self.allocator, cfg, req),
            .local => |cfg| localChat(self.allocator, cfg, req),
        };
    }
};

fn writeJSONEscaped(w: anytype, s: []const u8) !void {
    try w.writeAll("\"");
    for (s) |c| {
        switch (c) {
            '\\' => try w.writeAll("\\\\"),
            '"' => try w.writeAll("\\\""),
            else => try w.writeChar(c),
        }
    }
    try w.writeAll("\"");
}

fn buildChatBody(
    allocator: std.mem.Allocator,
    req: ChatRequest,
) ![]u8 {
    var buf = std.ArrayList(u8).init(allocator);
    const w = buf.writer();

    try w.writeAll("{");
    try w.writeAll("\"model\":");
    try writeJSONEscaped(w, req.model);
    try w.writeAll(",\"messages\":[");

    var first = true;
    for (req.messages) |msg| {
        if (!first) try w.writeAll(",");
        first = false;

        try w.writeAll("{\"role\":");
        try writeJSONEscaped(w, roleToStr(msg.role));
        try w.writeAll(",\"content\":");
        try writeJSONEscaped(w, msg.content);
        try w.writeAll("}");
    }

    try w.writeAll("]");
    if (req.temperature) |t| {
        try w.print(",\"temperature\":{d}", .{t});
    }
    try w.writeAll("}");

    return buf.toOwnedSlice();
}

fn roleToStr(r: Role) []const u8 {
    return switch (r) {
        .system => "system",
        .user => "user",
        .assistant => "assistant",
    };
}

fn doPost(
    allocator: std.mem.Allocator,
    url: []const u8,
    api_key: ?[]const u8,
    body: []const u8,
) ![]u8 {
    var client = std.http.Client{ .allocator = allocator };
    defer client.deinit();

    const uri = try std.Uri.parse(url);
    var req = try client.open(.POST, uri, .{
        .extra_headers = if (api_key) |key| &[_]std.http.Header{
            .{ .name = "Authorization", .value = "Bearer " ++ key },
            .{ .name = "Content-Type", .value = "application/json" },
        } else &[_]std.http.Header{
            .{ .name = "Content-Type", .value = "application/json" },
        },
    });
    defer req.deinit();

    try req.send(body);
    try req.finish();
    try req.wait();

    return try req.reader().readAllAlloc(
        allocator,
        1024 * 1024,
    );
}

fn parseOpenAIResponse(
    allocator: std.mem.Allocator,
    json_bytes: []const u8,
) !ChatResponse {
    var doc = try std.json.parseFromSlice(
        std.json.Value,
        allocator,
        json_bytes,
        .{},
    );
    defer doc.deinit();

    const root = doc.value;
    const content_val =
        root.object.get("choices").?
        .array.items[0]
        .object.get("message").?
        .object.get("content").?
        .string;

    return ChatResponse{ .text = content_val };
}

fn openaiChat(
    allocator: std.mem.Allocator,
    cfg: OpenAI,
    req: ChatRequest,
) !ChatResponse {
    const body = try buildChatBody(allocator, req);
    defer allocator.free(body);

    const raw = try doPost(
        allocator,
        cfg.base_url ++ "/v1/chat/completions",
        .{cfg.api_key},
        body,
    );

    return try parseOpenAIResponse(allocator, raw);
}

fn localChat(
    allocator: std.mem.Allocator,
    cfg: Local,
    req: ChatRequest,
) !ChatResponse {
    const body = try buildChatBody(allocator, req);
    defer allocator.free(body);

    // Ollama supports OpenAI compatible API under /v1/chat/completions
    const raw = try doPost(
        allocator,
        cfg.base_url ++ "/v1/chat/completions",
        null,
        body,
    );
    return try parseOpenAIResponse(allocator, raw);
}
