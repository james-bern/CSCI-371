template <typename T> struct StretchyBuffer {
    int length;
    int capacity;
    T *data;
        T &operator [](int index) { return data[index]; }
};

template <typename T> void sbuff_push_back(StretchyBuffer<T> *buffer, T element) {
    if (buffer->capacity == 0) {
        buffer->capacity = 16;
        buffer->data = (T *) malloc(buffer->capacity * sizeof(T));
    }
    if (buffer->length == buffer->capacity) {
        buffer->capacity *= 2;
        buffer->data = (T *) realloc(buffer->data, buffer->capacity * sizeof(T));
    }
    buffer->data[buffer->length++] = element;
}

template <typename T> void sbuff_free(StretchyBuffer<T> *buffer) {
    buffer->length = 0;
    buffer->capacity = 0;
    free(buffer->data);
    buffer->data = NULL;
}
